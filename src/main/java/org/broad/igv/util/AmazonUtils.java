package org.broad.igv.util;

import com.amazonaws.AmazonServiceException;
import com.amazonaws.HttpMethod;
import com.amazonaws.SdkClientException;
import com.amazonaws.auth.AWSStaticCredentialsProvider;
import com.amazonaws.auth.AnonymousAWSCredentials;
import com.amazonaws.auth.BasicSessionCredentials;
import com.amazonaws.regions.Regions;
import com.amazonaws.services.cognitoidentity.AmazonCognitoIdentity;
import com.amazonaws.services.cognitoidentity.AmazonCognitoIdentityClientBuilder;
import com.amazonaws.services.cognitoidentity.model.*;

import com.amazonaws.services.s3.AmazonS3;
import com.amazonaws.services.s3.AmazonS3ClientBuilder;
import com.amazonaws.services.s3.AmazonS3URI;
import com.amazonaws.services.s3.model.*;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import htsjdk.samtools.util.Tuple;
import org.apache.log4j.Logger;
import org.broad.igv.DirectoryManager;
import org.broad.igv.aws.S3Object;
import org.broad.igv.ga4gh.OAuthUtils;

import java.io.IOException;
import java.net.URL;
import java.util.ArrayList;
import java.util.Date;

public class AmazonUtils {
    private static Logger log = Logger.getLogger(AmazonUtils.class);

    // AWS specific objects
    private static AmazonS3 s3Client;
    private static JsonObject igv_oauth_conf = GetCognitoConfig();
    private static Credentials credentials;

    public static JsonObject GetCognitoConfig() {
        try {
            // Get AWS-specific Cognito details from on-disk JSON preferences
            String oauthConfig = DirectoryManager.getIgvDirectory() + "/oauth-config.json";
            JsonParser parser = new JsonParser();
            String json_file = FileUtils.getContents(oauthConfig);
            JsonObject json_obj = parser.parse(json_file).getAsJsonObject();
            return json_obj;
        } catch (IOException io) {
            log.error("Cannot read oauth-config.json config file");
        }

        return null;
    }

    /**
     * Returns the AWS credentials
     *
     * @return returns the credentials based on the access token returned from the user pool.
     */
    public static Credentials GetCognitoAWSCredentials() {

        JsonObject igv_oauth_conf = GetCognitoConfig();
        JsonObject response = OAuthUtils.getResponse();

        JsonObject payload = JWTParser.getPayload(response.get("id_token").getAsString());

        log.debug("JWT payload id token: "+payload);

        String idTokenStr = response.get("id_token").getAsString();
        String idProvider = payload.get("iss").toString().replace("https://", "").replace("\"", "");

        AnonymousAWSCredentials awsCreds = new AnonymousAWSCredentials();
        AmazonCognitoIdentity provider = AmazonCognitoIdentityClientBuilder
                .standard()
                .withCredentials(new AWSStaticCredentialsProvider(awsCreds))
                .withRegion(Regions.fromName(igv_oauth_conf.get("aws_region").getAsString()))
                .build();

        GetIdRequest idrequest = new GetIdRequest();
        idrequest.setIdentityPoolId(igv_oauth_conf.get("aws_cognito_fed_pool_id").getAsString());
        idrequest.addLoginsEntry(idProvider, idTokenStr);
        GetIdResult idResult = provider.getId(idrequest);

        GetCredentialsForIdentityRequest request = new GetCredentialsForIdentityRequest();
        request.setIdentityId(idResult.getIdentityId());
        request.addLoginsEntry(idProvider, idTokenStr);

        GetCredentialsForIdentityResult result = provider.getCredentialsForIdentity(request);

        credentials = result.getCredentials();

        return credentials;
    }


    /**
     * Makes sure the S3 client is available for bucket operations and/or generation of pre-signed urls
     *
     * @param credentials AWS credentials
     *
     */
    public static void updateS3Client(Credentials credentials) {
        BasicSessionCredentials awsCreds = new BasicSessionCredentials(credentials.getAccessKeyId(),
                credentials.getSecretKey(),
                credentials.getSessionToken());

        s3Client = AmazonS3ClientBuilder.standard()
                .withCredentials(new AWSStaticCredentialsProvider(awsCreds))
                .withRegion(Regions.fromName(igv_oauth_conf.get("aws_region").getAsString()))
                .build();
    }


    /**
     * This method returns the details of the user and bucket lists.
     * @return bucket list
     */
    public static ArrayList<String> ListBucketsForUser() {
        ArrayList<String> bucketsList = new ArrayList<>();

        for (Bucket bucket : s3Client.listBuckets()) {
            bucketsList.add(bucket.getName());
        }
        return bucketsList;
    }

    /**
     * Returns a list of bucket objects given a specific path.
     *
     * Adopted from:
     * https://docs.aws.amazon.com/AmazonS3/latest/dev/ListingObjectKeysUsingJava.html
     *
     * @return List of bucket objects
     */

    public static ArrayList<S3Object> ListBucketObjects(String bucketName, String prefix) {
        ArrayList<S3Object> objects = new ArrayList<>();
        log.debug("Listing objects for bucketName: "+ bucketName);
        try {
            // https://docs.aws.amazon.com/AmazonS3/latest/dev/UsingMetadata.html
            // """
            // The Amazon S3 data model is a flat structure: you create a bucket, and the bucket stores
            // objects. There is no hierarchy of subbuckets or subfolders; however, you can infer logical
            // hierarchy using key name prefixes and delimiters as the Amazon S3 console does. The Amazon
            // S3 console supports a concept of folders.
            // """
            ListObjectsV2Request req = new ListObjectsV2Request().withBucketName(bucketName).withPrefix(prefix).withDelimiter("/");
            ListObjectsV2Result result;

            do {
                result = s3Client.listObjectsV2(req);
                log.debug("S3 bucket prefix: "+result.getPrefix());

                for (String folder : result.getCommonPrefixes()) {
                    log.debug("S3 Bucket folder: "+folder);
                    folder = folder.substring(0, folder.length()-1); // Chop off last / of the folder for UI purposes
                    objects.add(new S3Object(folder.replace(prefix, ""), true));
                }

                for (S3ObjectSummary objectSummary : result.getObjectSummaries()) {
                    log.debug("S3 Bucket key: "+objectSummary.getKey());
                    objects.add(new S3Object(objectSummary.getKey().replace(prefix, ""), false));
                }
                // If there are more than maxKeys keys in the bucket, get a continuation token
                // and list the next objects.
                String token = result.getNextContinuationToken();
                log.debug("Next S3 bucket pagination continuation Token: " + token);
                req.setContinuationToken(token);
            } while (result.isTruncated());

        } catch(AmazonServiceException e) {
            // The call was transmitted successfully, but Amazon S3 couldn't process
            // it, so it returned an error response.
            e.printStackTrace();
        } catch(SdkClientException e) {
            // Amazon S3 couldn't be contacted for a response, or the client
            // couldn't parse the response from Amazon S3.
            e.printStackTrace();
        }

        return objects;
    }

    public static Tuple<String, String> bucketAndKey(String S3urlString) {
        AmazonS3URI s3URI = new AmazonS3URI(S3urlString);

        return new Tuple(s3URI.getBucket(), s3URI.getKey());
    }

    /**
     * Generates a so-called "pre-signed" URLs (https://docs.aws.amazon.com/AmazonS3/latest/dev/ShareObjectPreSignedURLJavaSDK.html)
     * such as:
     *
     * s3://igv-bam-test/NA12878.bam
     * https://s3-<REGION>.amazonaws.com/igv-bam-test/NA12878.bam
     *
     * @param bucketName
     * @param objectKey
     * @return
     */
    public static URL translateAmazonCloudURL(String bucketName, String objectKey, Date expirationTime) {
        // We generate presigned URLs out of the S3 bucket because loadTracks->HttpUtils does not understand s3://
        // ... good old IETF rfc2396 enforcing standard protocol schemes.
        // This wrapper name is consistent with GoogleUtils similarly-named class method.

        log.debug("Generating pre-signed URL for: "+ bucketName + "/" + objectKey);

        GeneratePresignedUrlRequest generatePresignedUrlRequest =
                new GeneratePresignedUrlRequest(bucketName, objectKey)
                        .withMethod(HttpMethod.GET)
                        .withExpiration(expirationTime);

        return s3Client.generatePresignedUrl(generatePresignedUrlRequest);
    }

    public static String translateAmazonCloudURL(String urlString) {
        Tuple<String, String> bandk = bucketAndKey(urlString);
        String bucketName = bandk.a;
        String objectKey = bandk.b;

        log.debug("Generating pre-signed URL for: "+ bandk.a + "/" + bandk.b);

        GeneratePresignedUrlRequest generatePresignedUrlRequest =
                new GeneratePresignedUrlRequest(bucketName, objectKey)
                        .withMethod(HttpMethod.GET)
                        .withExpiration(OAuthUtils.getExpirationDate());

        return s3Client.generatePresignedUrl(generatePresignedUrlRequest).toString();
    }

}
