package org.broad.igv.util;

import com.google.gson.JsonObject;
import org.broad.igv.Globals;
import org.broad.igv.aws.IGVS3Object;
import org.broad.igv.google.OAuthProvider;
import org.broad.igv.google.OAuthUtils;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.ui.IGVMenuBar;
import software.amazon.awssdk.auth.credentials.AnonymousCredentialsProvider;
import software.amazon.awssdk.auth.credentials.AwsSessionCredentials;
import software.amazon.awssdk.auth.credentials.DefaultCredentialsProvider;
import software.amazon.awssdk.auth.credentials.StaticCredentialsProvider;
import software.amazon.awssdk.core.exception.SdkClientException;
import software.amazon.awssdk.core.exception.SdkServiceException;
import software.amazon.awssdk.regions.Region;
import software.amazon.awssdk.regions.providers.DefaultAwsRegionProviderChain;
import software.amazon.awssdk.services.cognitoidentity.CognitoIdentityClient;
import software.amazon.awssdk.services.cognitoidentity.CognitoIdentityClientBuilder;
import software.amazon.awssdk.services.cognitoidentity.model.GetIdRequest;
import software.amazon.awssdk.services.cognitoidentity.model.GetIdResponse;
import software.amazon.awssdk.services.cognitoidentity.model.GetOpenIdTokenRequest;
import software.amazon.awssdk.services.cognitoidentity.model.GetOpenIdTokenResponse;
import software.amazon.awssdk.services.s3.S3Client;
import software.amazon.awssdk.services.s3.model.*;
import software.amazon.awssdk.services.s3.paginators.ListObjectsV2Iterable;
import software.amazon.awssdk.services.sts.StsClient;
import software.amazon.awssdk.services.sts.model.AssumeRoleWithWebIdentityRequest;
import software.amazon.awssdk.services.sts.model.AssumeRoleWithWebIdentityResponse;
import software.amazon.awssdk.services.sts.model.Credentials;

import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.net.URI;
import java.net.URL;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.concurrent.CompletableFuture;
import java.util.stream.Collectors;

public class AmazonUtils {
    private static Logger log = LogManager.getLogger(AmazonUtils.class);

    // AWS specific objects
    private static S3Client s3Client;
    private static List<String> bucketsFinalList = new ArrayList<>();
    private static CognitoIdentityClient cognitoIdentityClient;
    private static Region AWSREGION;
    private static Boolean awsCredentialsPresent = null;


    /**
     * Maps s3:// URLs to presigned URLs
     */
    private static Map<String, String> s3ToPresignedMap = new HashMap<>();

    /**
     * Maps aws presigned URLs to s3://.  This is needed in some cases (e.g. Tribble) to regenerate an expired URL
     */
    private static Map<String, String> presignedToS3Map = new HashMap<>();

    private static JsonObject CognitoConfig;

    public static void setCognitoConfig(JsonObject json) {
        CognitoConfig = json;
        awsCredentialsPresent = CognitoConfig.get("auth_provider").getAsString().contains("Amazon");
        if (IGVMenuBar.getInstance() != null) {
            IGVMenuBar.getInstance().updateAWSMenu();
        }
    }

    public static JsonObject GetCognitoConfig() {
        return CognitoConfig;
    }


    /**
     * Test to see if aws credentials are avaialable, either through IGV configuration of Cognito, or from the
     * default provider chain.  See https://docs.aws.amazon.com/sdk-for-java/v1/developer-guide/credentials.html#credentials-default
     *
     * @return
     */
    public static boolean isAwsProviderPresent() {

        if (awsCredentialsPresent == null) {
            if (GetCognitoConfig() != null) {
                try {
                    if (GetCognitoConfig().get("auth_provider").getAsString().contains("Amazon")) {
                        log.info("AWS configuration found. AWS support enabled.");
                        awsCredentialsPresent = true;
                    } else {
                        log.info("AWS configuration not found.");
                        awsCredentialsPresent = false;
                    }
                } catch (NullPointerException np) {
                    awsCredentialsPresent = false;
                }
            } else {
                try {
                    DefaultCredentialsProvider.create().resolveCredentials();
                    log.info("AWS default credentials found. AWS support enabled.");
                    awsCredentialsPresent = true;
                } catch (Exception e) {
                    // This is an expected condition if credentials can't be found => don't log
                    awsCredentialsPresent = false;
                }
            }
        }
        return awsCredentialsPresent;
    }

    private static Region getAWSREGION() {

        if (AWSREGION == null) {
            if (GetCognitoConfig() != null) {
                AWSREGION = Region.of(GetCognitoConfig().get("aws_region").getAsString());
            } else {
                // TODO -- find region in default place
                try {
                    AWSREGION = (new DefaultAwsRegionProviderChain()).getRegion();
                } catch (Exception e) {
                    AWSREGION = Region.US_EAST_1;
                }
            }
        }
        return AWSREGION;
    }

    /**
     * Returns the AWS credentials
     *
     * @return returns the credentials based on the AWS STS access token returned from the AWS Cognito user pool.
     */
    public static Credentials GetCognitoAWSCredentials() {

        OAuthProvider provider = OAuthUtils.getInstance().getProvider("Amazon");

        JsonObject igv_oauth_conf = GetCognitoConfig();
        JsonObject response = provider.getResponse();

        // Handle non-user initiated S3 auth (IGV early startup), i.e user-specified GenomesLoader
        if (response == null) {
            // Go back to auth flow, not auth'd yet
            checkLogin();
            response = provider.getResponse();
        }

        JsonObject payload = JWTParser.getPayload(response.get("id_token").getAsString());

        log.debug("JWT payload id token: " + payload);

        // Collect necessary information from federated IdP for Authentication purposes
        String idTokenStr = response.get("id_token").getAsString();
        String idProvider = payload.get("iss").toString().replace("https://", "")
                .replace("\"", "");
        String email = payload.get("email").getAsString();
        String federatedPoolId = igv_oauth_conf.get("aws_cognito_fed_pool_id").getAsString();
        String cognitoRoleARN = igv_oauth_conf.get("aws_cognito_role_arn").getAsString();

        HashMap<String, String> logins = new HashMap<>();
        logins.put(idProvider, idTokenStr);

        // Avoid "software.amazon.awssdk.core.exception.SdkClientException: Unable to load credentials from any of the providers in the chain AwsCredentialsProviderChain("
        // The use of the AnonymousCredentialsProvider essentially bypasses the provider chain's requirement to access ~/.aws/credentials.
        // https://stackoverflow.com/questions/36604024/sts-saml-and-java-sdk-unable-to-load-aws-credentials-from-any-provider-in-the-c
        AnonymousCredentialsProvider anoCredProv = AnonymousCredentialsProvider.create();

        // https://sdk.amazonaws.com/java/api/latest/software/amazon/awssdk/services/cognitoidentity/CognitoIdentityClient.html
        // Build the Cognito client
        CognitoIdentityClientBuilder cognitoIdentityBuilder = CognitoIdentityClient.builder();

        cognitoIdentityBuilder.region(getAWSREGION()).credentialsProvider(anoCredProv);
        cognitoIdentityClient = cognitoIdentityBuilder.build();


        // https://docs.aws.amazon.com/cognito/latest/developerguide/authentication-flow.html
        // Basic (Classic) Authflow
        // Uses AssumeRoleWithWebIdentity and facilitates CloudTrail org.broad.igv.logging. Uses one more request but provides user traceability.
        GetIdRequest.Builder idrequest = GetIdRequest.builder().identityPoolId(federatedPoolId)
                .logins(logins);
        GetIdResponse idResult = cognitoIdentityClient.getId(idrequest.build());

        GetOpenIdTokenRequest.Builder openidrequest = GetOpenIdTokenRequest.builder().logins(logins).identityId(idResult.identityId());
        GetOpenIdTokenResponse openId = cognitoIdentityClient.getOpenIdToken(openidrequest.build());


        AssumeRoleWithWebIdentityRequest.Builder webidrequest = AssumeRoleWithWebIdentityRequest.builder().webIdentityToken(openId.token())
                .roleSessionName(email)
                .roleArn(cognitoRoleARN);

        AssumeRoleWithWebIdentityResponse stsClientResponse = StsClient.builder().credentialsProvider(anoCredProv)
                .region(getAWSREGION())
                .build()
                .assumeRoleWithWebIdentity(webidrequest.build());

//      // Enhanced (Simplified) Authflow
//      // Major drawback: Does not store federated user information on CloudTrail only authenticated role name appears in logs.
//
//        // "To provide end-user credentials, first make an unsigned call to GetId."
//        GetIdRequest.Builder idrequest = GetIdRequest.builder().identityPoolId(federatedPoolId)
//                                                               .logins(logins);
//        GetIdResponse idResult = cognitoIdentityClient.getId(idrequest.build());
//
//        // "Next, make an unsigned call to GetCredentialsForIdentity."
//        GetCredentialsForIdentityRequest.Builder authedIds = GetCredentialsForIdentityRequest.builder();
//        authedIds.identityId(idResult.identityId()).logins(logins);
//
//        GetCredentialsForIdentityResponse authedRes = cognitoIdentityClient.getCredentialsForIdentity(authedIds.build());
//
//        return authedRes.credentials()

        return stsClientResponse.credentials();
    }


    /**
     * Makes sure the S3 client is available for bucket operations and/or generation of pre-signed urls
     *
     * @param credentials AWS credentials
     */
    public static void updateS3Client(Credentials credentials) {
        final Region region = getAWSREGION();
        if (credentials == null) {
            s3Client = S3Client.builder().region(region).build();
        } else {
            AwsSessionCredentials creds = AwsSessionCredentials.create(credentials.accessKeyId(),
                    credentials.secretAccessKey(),
                    credentials.sessionToken());

            StaticCredentialsProvider s3CredsProvider = StaticCredentialsProvider.create(creds);

            s3Client = S3Client.builder().credentialsProvider(s3CredsProvider).region(region).build();
        }
    }


    /**
     * This method returns the details of the user and bucket lists.
     *
     * @return bucket list
     */
    public static List<String> ListBucketsForUser() {
        if (bucketsFinalList.isEmpty()) {

            if (GetCognitoConfig() != null) {
                OAuthUtils.getInstance().getProvider("Amazon").getAccessToken();
                updateS3Client(GetCognitoAWSCredentials());
            } else {
                updateS3Client(null);
            }
            List<String> bucketsList = new ArrayList<>();

            ListBucketsRequest listBucketsRequest = ListBucketsRequest.builder().build();
            ListBucketsResponse listBucketsResponse = s3Client.listBuckets(listBucketsRequest);
            listBucketsResponse.buckets().stream().forEach(x -> bucketsList.add(x.name()));

            // Filter out buckets that the user does not have permissions for
            bucketsFinalList = getReadableBuckets(bucketsList);

        }
        return bucketsFinalList;
    }

    public static HeadObjectResponse getObjectMetadata(String bucket, String key) {
        HeadObjectRequest HeadObjReq = HeadObjectRequest.builder()
                .bucket(bucket)
                .key(key).build();
        HeadObjectResponse HeadObjRes = s3Client.headObject(HeadObjReq);
        log.debug("getObjectMetadata(): "+HeadObjRes.toString());
        return HeadObjRes;
    }


    // Holds whether a S3 object is accessible or not and reason/error msg in case it's not.
    public static class s3ObjectAccessResult {
        private boolean objAvailable;
        private String errorReason;

        public boolean isObjectAvailable() {
            return objAvailable;
        }

        public void setObjAvailable(boolean objAvailable) {
            this.objAvailable = objAvailable;
        }

        public String getErrorReason() {
            return errorReason;
        }

        public void setErrorReason(String errorReason) {
            this.errorReason = errorReason;
        }
    }

    // Determines whether the object is immediately available.
    // On AWS this means present in STANDARD, STANDARD_IA, INTELLIGENT_TIERING object access tiers.
    // Tiers GLACIER and DEEP_ARCHIVE are not immediately retrievable without action.
    public static s3ObjectAccessResult isObjectAccessible(String bucket, String key) {
        // TODO:
        //  1. Determine if it's a public S3 resource first of all? If so, none of the logic below is needed
        s3ObjectAccessResult res = new s3ObjectAccessResult();

        // Safeguard for null corner case(s), assume we can access the object
        res.setObjAvailable(true);
        //res.setErrorReason("Object not found, perhaps a new tier was introduced on AWS?"); // not really an error

        HeadObjectResponse s3Meta;      // Head metadata from the S3 object
        String s3ObjectStorageStatus;   // Can it be retrieved immediately or not?
        String s3ObjectStorageClass;    // Which AWS S3 tier is this object in?

        // Simple "null" case. The object is directly accessible in
        // STANDARD, INFREQUENT_ACCESS, INTELLIGENT_TIERING
        // or any other "immediately available" tier.
        s3Meta = AmazonUtils.getObjectMetadata(bucket, key);
        if (s3Meta.storageClass() == null) {
            res.setErrorReason("Object is in an accessible tier, no errors are expected");
            res.setObjAvailable(true);
            return res; // nothing else to check, return early
        }

        // Determine in which state this object really is:
        // 1. Archived.
        // 2. In the process of being restored.
        // 3. Restored
        //
        // This is important because after restoration the object mantains the Tier (DEEP_ARCHIVE) instead of
        // transitioning that attribute to STANDARD, we must look at head_object response for the "Restore"
        // attribute.
        //
        // Possible error reason messages for the users are:

        s3ObjectStorageClass = s3Meta.storageClass().toString();

        String archived = "Amazon S3 object is in " + s3ObjectStorageClass + " storage tier, not accessible at this moment. " +
                "Please contact your local system administrator about object: s3://" + bucket + "/" + key;
        String restoreInProgress = "Amazon S3 object is in " + s3ObjectStorageClass + " and being restored right now, please be patient, this can take up to 48h. " +
                "For further enquiries about this dataset, please use the following path when communicating with your system administrator: s3://" + bucket + "/" + key;

        if (s3ObjectStorageClass.contains("DEEP_ARCHIVE") ||
                s3ObjectStorageClass.contains("GLACIER")) {
            try {
                s3ObjectStorageStatus = s3Meta.sdkHttpResponse().headers().get("x-amz-restore").toString();
            } catch (NullPointerException npe) {
                res.setObjAvailable(false);
                res.setErrorReason(archived);
                return res;
            }

            if (s3ObjectStorageStatus.contains("ongoing-request=\"true\"")) {
                res.setObjAvailable(false);
                res.setErrorReason(restoreInProgress);

                // "If an archive copy is already restored, the header value indicates when Amazon S3 is scheduled to delete the object copy"
            } else if (s3ObjectStorageStatus.contains("ongoing-request=\"false\"") && s3ObjectStorageStatus.contains("expiry-date=")) {
                res.setObjAvailable(true);
            } else {
                // The object has never been restored?
                res.setObjAvailable(false);
                res.setErrorReason(archived);
            }
        }

        return res;
    }

    private static List<String> getReadableBuckets(List<String> buckets) {
        List<CompletableFuture<String>> futures =
                buckets.stream()
                        .map(bucket -> {
                            CompletableFuture<String> future = CompletableFuture.supplyAsync(() -> {
                                if (AmazonUtils.ListBucketObjects(bucket, "").size() > 0) {
                                    return bucket;
                                }

                                return null;
                            });

                            return future;
                        })
                        .collect(Collectors.toList());

        List<String> result =
                futures.stream()
                        .map(CompletableFuture::join)
                        .collect(Collectors.toList());

        result.removeAll(Collections.singleton(null));

        return result;
    }

    /**
     * Returns a list of bucket objects given a specific path.
     * <p>
     * Adopted from:
     * https://docs.aws.amazon.com/AmazonS3/latest/dev/ListingObjectKeysUsingJava.html
     *
     * @return List of bucket objects
     */

    public static ArrayList<IGVS3Object> ListBucketObjects(String bucketName, String prefix) {
        ArrayList<IGVS3Object> objects = new ArrayList<>();
        log.debug("Listing objects for bucketName: " + bucketName);

        if (GetCognitoConfig() != null) {
            OAuthUtils.getInstance().getProvider("Amazon").getAccessToken();
            updateS3Client(GetCognitoAWSCredentials());
        } else {
            updateS3Client(null);
        }

        try {
            // https://docs.aws.amazon.com/AmazonS3/latest/dev/UsingMetadata.html
            // """
            // The Amazon S3 data model is a flat structure: you create a bucket, and the bucket stores
            // objects. There is no hierarchy of subbuckets or subfolders; however, you can infer logical
            // hierarchy using key name prefixes and delimiters as the Amazon S3 console does. The Amazon
            // S3 console supports a concept of folders.
            // """
            ListObjectsV2Request listReq = ListObjectsV2Request.builder().bucket(bucketName)
                    .prefix(prefix)
                    .delimiter("/")
                    .build();

            ListObjectsV2Response response = s3Client.listObjectsV2(listReq);
            ListObjectsV2Iterable resultIt = s3Client.listObjectsV2Paginator(listReq);
            String folder_prefix;

            do {
                for (CommonPrefix folder : resultIt.commonPrefixes()) {
                    log.debug("S3 Bucket folder: " + folder);
                    folder_prefix = folder.prefix().substring(0, folder.prefix().length() - 1); // Chop off last / of the folder for UI purposes
                    objects.add(new IGVS3Object(folder_prefix.replace(prefix, ""), true, "STANDARD"));
                }

                for (S3Object content : response.contents()) {
                    log.debug("S3 Bucket key: " + content.key());
                    objects.add(new IGVS3Object(content.key().replace(prefix, ""), false, content.storageClassAsString()));
                }
                // If there are more than maxKeys keys in the bucket, get a continuation token
                // and list the next objects.
                String token = response.nextContinuationToken();
                log.debug("Next S3 bucket pagination continuation Token: " + token);
                response.continuationToken();
            } while (response.isTruncated());

        } catch (SdkServiceException | SdkClientException e) {
            // The call was transmitted successfully, but Amazon S3 couldn't process
            // it, so it returned an error response.
            log.debug("AccessDenied for ListBucket " + bucketName + " with prefix " + prefix);
        }

        return objects;
    }

    public static String getBucketFromS3URL(String s3URL) {
        AmazonS3URI s3URI = new AmazonS3URI(s3URL);
        return s3URI.getBucket();

    }

    public static String getKeyFromS3URL(String s3URL) {
        AmazonS3URI s3URI = new AmazonS3URI(s3URL);
        return s3URI.getKey();
    }

    // Amazon S3 Presign URLs
    // Also keeps an internal mapping between ResourceLocator and active/valid signed URLs.
    private static String createPresignedURL(String s3Path) throws IOException {
        // TODO: Ideally the presigned URL should be generated without any of the Cognito being involved first?
        // Make sure access token are valid (refreshes token internally)

        S3Presigner s3Presigner;

        if (GetCognitoConfig() != null) {
            OAuthProvider provider = OAuthUtils.getInstance().getProvider("Amazon");
            provider.getAccessToken();

            Credentials credentials = GetCognitoAWSCredentials();
            AwsSessionCredentials creds = AwsSessionCredentials.create(credentials.accessKeyId(),
                    credentials.secretAccessKey(),
                    credentials.sessionToken());
            StaticCredentialsProvider awsCredsProvider = StaticCredentialsProvider.create(creds);

            s3Presigner = S3Presigner.builder()
                    .expiration(provider.getExpirationTime())       // Duration.ofSeconds(30)  // <= for testing
                    .awsCredentials(awsCredsProvider)
                    .region(getAWSREGION())
                    .build();
        } else {
            s3Presigner = S3Presigner.builder().build();
        }

        String bucket = getBucketFromS3URL(s3Path);
        String key = getKeyFromS3URL(s3Path);

        URI presigned = s3Presigner.presignS3DownloadLink(bucket, key);
        log.debug("AWS presigned URL from translateAmazonCloudURL is: " + presigned);
        return presigned.toString();
    }

    /**
     * @param s3UrlString
     * @return
     * @throws IOException
     */

    public static String translateAmazonCloudURL(String s3UrlString) throws IOException {
        String presignedUrl = s3ToPresignedMap.get(s3UrlString);
        if (presignedUrl == null || !isPresignedURLValid(new URL(presignedUrl))) {
            presignedUrl = createPresignedURL(s3UrlString);
            s3ToPresignedMap.put(s3UrlString, presignedUrl);
            presignedToS3Map.put(presignedUrl, s3UrlString);
        }
        return presignedUrl;
    }

    public static Boolean isAwsS3Path(String path) {
        // TODO: perhaps add some more checks
        return (path.startsWith("s3://"));
    }

    public static boolean isPresignedURL(String urlString) {
        return presignedToS3Map.containsKey(urlString);
    }

    public static String updatePresignedURL(String urlString) throws IOException {
        String s3UrlString = presignedToS3Map.get(urlString);
        if(s3UrlString == null) {
            throw new RuntimeException("Unrecognized presigned url: " + urlString);
        } else {
            return translateAmazonCloudURL(s3UrlString);
        }
    }

    /**
     * If using Cognito, check that the use is logged in, and prompt for login if not.
     */
    public static void checkLogin() {
        if (GetCognitoConfig() != null &&
                !OAuthUtils.getInstance().getProvider("Amazon").isLoggedIn()) {
            OAuthUtils.getInstance().getProvider("Amazon").doSecureLogin();
        }
    }

    /**
     * Checks whether a (pre)signed url is still accessible or it has expired, offline.
     * No extra request/head is required to the presigned object since we have all information
     * available on the AWS URL parameters themselves, namely:
     * <p>
     * X-Amz-Expires=12 (in seconds)
     * X-Amz-Date=20190725T045535Z
     * <p>
     * NOTE: X-Amz-Date is expressed in Zulu (military) time. The rest is on UTC, so we'll use UTC
     **/

    private static boolean isPresignedURLValid(URL url) {
        boolean isValidSignedUrl;

        try {
            long presignedTime = signedURLValidity(url);
            isValidSignedUrl = presignedTime - System.currentTimeMillis() - Globals.TOKEN_EXPIRE_GRACE_TIME > 0; // Duration in milliseconds
            if(!isValidSignedUrl) {
                System.out.println("URL expired: " + url.toExternalForm());
            }
        } catch (ParseException e) {
            log.error("The AWS signed URL date parameter X-Amz-Date has incorrect formatting");
            isValidSignedUrl = false;
        } catch (UnsupportedEncodingException e) {
            e.printStackTrace();
            isValidSignedUrl = false;
        }

        return isValidSignedUrl;
    }

    private static long signedURLValidity(URL url) throws ParseException, UnsupportedEncodingException {
        Map<String, String> params = StringUtils.splitQuery(url);
        String amzDateStr = params.get("X-Amz-Date");
        long amzExpires = Long.parseLong(params.get("X-Amz-Expires"));

        SimpleDateFormat formatter = new SimpleDateFormat("yyyyMMdd'T'HHmmss'Z'");
        formatter.setTimeZone(TimeZone.getTimeZone("UTC")); // Z(ulu) -> UTC
        Date amzDate = formatter.parse(amzDateStr);

        long timeOfExpirationMillis = amzDate.getTime() + amzExpires * 1000;

        log.debug("The date of expiration is " + amzDate + ", expires after " + amzExpires + " seconds for url: " + url);
        return timeOfExpirationMillis;
    }
}
