package org.broad.igv.util;

import com.amazonaws.auth.AWSStaticCredentialsProvider;
import com.amazonaws.auth.AnonymousAWSCredentials;
import com.amazonaws.auth.BasicSessionCredentials;
import com.amazonaws.regions.Regions;
import com.amazonaws.services.cognitoidentity.AmazonCognitoIdentity;
import com.amazonaws.services.cognitoidentity.AmazonCognitoIdentityClientBuilder;
import com.amazonaws.services.cognitoidentity.model.*;

import com.amazonaws.services.s3.AmazonS3;
import com.amazonaws.services.s3.AmazonS3ClientBuilder;
import com.amazonaws.services.s3.model.Bucket;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import org.apache.log4j.Logger;
import org.broad.igv.DirectoryManager;

import java.io.IOException;

import static com.amazonaws.auth.profile.internal.ProfileKeyConstants.REGION;

public class AmazonUtils {
    private static Logger log = Logger.getLogger(AmazonUtils.class);

    private static JsonObject GetCognitoConfig() {
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
     * @param response         to get the username for the login map.
     * @return returns the credentials based on the access token returned from the user pool.
     */
    public static Credentials GetCognitoAWSCredentials(JsonObject response) {

        String response_str = response.toString();

        JsonObject igv_oauth_conf = GetCognitoConfig();
        JsonObject payload = JWTParser.getPayload(response.get("id_token").getAsString());
        JsonObject payload_access = JWTParser.getPayload(response.get("access_token").getAsString());

        log.info("JWT payload id token: "+payload);
        log.info("JWT payload access token: "+payload_access);

        String id = response_str;
        String idprovider = payload.get("iss").toString().replace("https://", "");

        AnonymousAWSCredentials awsCreds = new AnonymousAWSCredentials();
        AmazonCognitoIdentity provider = AmazonCognitoIdentityClientBuilder
                .standard()
                .withCredentials(new AWSStaticCredentialsProvider(awsCreds))
                .withRegion(Regions.fromName(igv_oauth_conf.get("aws_region").getAsString()))
                .build();

        GetIdRequest idrequest = new GetIdRequest();
        idrequest.setIdentityPoolId(igv_oauth_conf.get("aws_cognito_fed_pool_id").getAsString());
        idrequest.addLoginsEntry(idprovider, id);
        GetIdResult idResult = provider.getId(idrequest);

        GetCredentialsForIdentityRequest request = new GetCredentialsForIdentityRequest();
        request.setIdentityId(idResult.getIdentityId());
        request.addLoginsEntry(idprovider, id);

        GetCredentialsForIdentityResult result = provider.getCredentialsForIdentity(request);
        return result.getCredentials();
    }


     /*
     * This method returns the details of the user and bucket lists.
     *
     * @param credentials Credentials to be used for displaying buckets
     * @return
     */
    public static String ListBucketsForUser(Credentials credentials) {
        JsonObject igv_oauth_conf = GetCognitoConfig();

        BasicSessionCredentials awsCreds = new BasicSessionCredentials(credentials.getAccessKeyId(),
                                                                       credentials.getSecretKey(),
                                                                       credentials.getSessionToken());
        AmazonS3 s3Client = AmazonS3ClientBuilder.standard()
                .withCredentials(new AWSStaticCredentialsProvider(awsCreds))
                .withRegion(Regions.fromName(igv_oauth_conf.get("aws_region").getAsString()))
                .build();
        StringBuilder bucketslist = new StringBuilder();

        bucketslist.append("===========Credentials Details.=========== \n");
        bucketslist.append("Accesskey = " + credentials.getAccessKeyId() + "\n");
        bucketslist.append("Secret = " + credentials.getSecretKey() + "\n");
        bucketslist.append("SessionToken = " + credentials.getSessionToken() + "\n");
        bucketslist.append("============Bucket Lists===========\n");

        for (Bucket bucket : s3Client.listBuckets()) {
            bucketslist.append(bucket.getName());
            bucketslist.append("\n");

            System.out.println(" - " + bucket.getName());
        }
        return bucketslist.toString();
    }

}
