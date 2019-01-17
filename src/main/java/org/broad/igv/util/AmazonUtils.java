package org.broad.igv.util;

import com.amazonaws.HttpMethod;
import com.amazonaws.auth.AWSStaticCredentialsProvider;
import com.amazonaws.auth.AnonymousAWSCredentials;
import com.amazonaws.auth.BasicSessionCredentials;
import com.amazonaws.regions.Regions;
import com.amazonaws.services.cognitoidentity.AmazonCognitoIdentity;
import com.amazonaws.services.cognitoidentity.AmazonCognitoIdentityClientBuilder;
import com.amazonaws.services.cognitoidentity.model.*;

import com.amazonaws.services.s3.AmazonS3;
import com.amazonaws.services.s3.AmazonS3Client;
import com.amazonaws.services.s3.AmazonS3ClientBuilder;
import com.amazonaws.services.s3.model.Bucket;
import com.amazonaws.services.s3.model.GeneratePresignedUrlRequest;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import org.apache.log4j.Logger;
import org.broad.igv.DirectoryManager;

import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.Date;

public class AmazonUtils {
    private static Logger log = Logger.getLogger(AmazonUtils.class);
    private static AmazonS3 s3Client;



    // XXX: A similar method must exist for this, use that one instead and/or move to more general class
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
     * @param response contains all the OAuth/OIDC information required to generate AWS credentials
     * @return returns the credentials based on the access token returned from the user pool.
     */
    public static Credentials GetCognitoAWSCredentials(JsonObject response) {

        JsonObject igv_oauth_conf = GetCognitoConfig();
        JsonObject payload = JWTParser.getPayload(response.get("id_token").getAsString());

        log.debug("JWT payload id token: "+payload);

        String id_token_str = response.get("id_token").getAsString();
        String idprovider = payload.get("iss").toString().replace("https://", "").replace("\"", "");

        AnonymousAWSCredentials awsCreds = new AnonymousAWSCredentials();
        AmazonCognitoIdentity provider = AmazonCognitoIdentityClientBuilder
                .standard()
                .withCredentials(new AWSStaticCredentialsProvider(awsCreds))
                .withRegion(Regions.fromName(igv_oauth_conf.get("aws_region").getAsString()))
                .build();

        GetIdRequest idrequest = new GetIdRequest();
        idrequest.setIdentityPoolId(igv_oauth_conf.get("aws_cognito_fed_pool_id").getAsString());
        idrequest.addLoginsEntry(idprovider, id_token_str);
        GetIdResult idResult = provider.getId(idrequest);

        GetCredentialsForIdentityRequest request = new GetCredentialsForIdentityRequest();
        request.setIdentityId(idResult.getIdentityId());
        request.addLoginsEntry(idprovider, id_token_str);

        GetCredentialsForIdentityResult result = provider.getCredentialsForIdentity(request);
        return result.getCredentials();
    }


     /*
     * This method returns the details of the user and bucket lists.
     *
     * @param credentials Credentials to be used for displaying buckets
     * @return
     */
    public String ListBucketsForUser(Credentials credentials) {
        // XXX: Fix the region and oauth-config.json stuff better
        JsonObject igv_oauth_conf = GetCognitoConfig();

        // XXX: Is converting cognito credentials to plain AWS creds really needed?
        BasicSessionCredentials awsCreds = new BasicSessionCredentials(credentials.getAccessKeyId(),
                                                                       credentials.getSecretKey(),
                                                                       credentials.getSessionToken());
        s3Client = AmazonS3ClientBuilder.standard()
                .withCredentials(new AWSStaticCredentialsProvider(awsCreds))
                .withRegion(Regions.fromName(igv_oauth_conf.get("aws_region").getAsString()))
                .build();

        StringBuilder bucketslist = new StringBuilder();

        for (Bucket bucket : s3Client.listBuckets()) {
            bucketslist.append(bucket.getName());
            bucketslist.append("\n");

            System.out.println(" - " + bucket.getName());
        }
        return bucketslist.toString();
    }


    public static void updateS3Client(Credentials credentials) {
        // XXX: Fix the region and oauth-config.json stuff better
        JsonObject igv_oauth_conf = GetCognitoConfig();

        // XXX: Is converting cognito credentials to plain AWS creds really needed?
        BasicSessionCredentials awsCreds = new BasicSessionCredentials(credentials.getAccessKeyId(),
                credentials.getSecretKey(),
                credentials.getSessionToken());

        s3Client = AmazonS3ClientBuilder.standard()
                .withCredentials(new AWSStaticCredentialsProvider(awsCreds))
                .withRegion(Regions.fromName(igv_oauth_conf.get("aws_region").getAsString()))
                .build();
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
        // We generate presigned URLs out of the S3 bucket because loadTracks->HttpUtils do not understand s3://
        // ... good old IETF rfc2396 enforcing standard protocol schemes.
        // This wrapper name is consistent with GoogleUtils similarly-named class method.

        // XXX: Make sure expiry/refresh times do not affect pre-signed URLs
        // XXX: equivalent to os.path.sep in Java?
        log.debug("Generating pre-signed URL for: "+ bucketName + "/" + objectKey);

        GeneratePresignedUrlRequest generatePresignedUrlRequest =
                new GeneratePresignedUrlRequest(bucketName, objectKey)
                        .withMethod(HttpMethod.GET)
                        .withExpiration(expirationTime);

        return s3Client.generatePresignedUrl(generatePresignedUrlRequest);
    }
}
