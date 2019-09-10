package org.broad.igv.util;

/*
*
* This is a transitional class until the official java-aws-sdk-v2 includes a S3 URL presigners class, see:
*
* https://github.com/aws/aws-sdk-java-v2/issues/849#issuecomment-468892839
* https://github.com/aws/aws-sdk-java-v2/issues/203
*
*/

import java.io.ByteArrayInputStream;
import java.net.URI;
import java.time.Duration;
import java.time.Instant;

import software.amazon.awssdk.auth.credentials.AwsCredentialsProvider;
import software.amazon.awssdk.auth.credentials.DefaultCredentialsProvider;
import software.amazon.awssdk.auth.signer.AwsS3V4Signer;
import software.amazon.awssdk.auth.signer.AwsSignerExecutionAttribute;
import software.amazon.awssdk.core.ResponseInputStream;
import software.amazon.awssdk.core.client.config.ClientOverrideConfiguration;
import software.amazon.awssdk.core.exception.SdkClientException;
import software.amazon.awssdk.core.interceptor.Context.BeforeTransmission;
import software.amazon.awssdk.core.interceptor.ExecutionAttributes;
import software.amazon.awssdk.core.interceptor.ExecutionInterceptor;
import software.amazon.awssdk.core.sync.RequestBody;
import software.amazon.awssdk.http.AbortableInputStream;
import software.amazon.awssdk.http.ExecutableHttpRequest;
import software.amazon.awssdk.http.HttpExecuteRequest;
import software.amazon.awssdk.http.HttpExecuteResponse;
import software.amazon.awssdk.http.SdkHttpClient;
import software.amazon.awssdk.http.SdkHttpFullRequest;
import software.amazon.awssdk.http.SdkHttpResponse;
import software.amazon.awssdk.regions.Region;
import software.amazon.awssdk.regions.providers.DefaultAwsRegionProviderChain;
import software.amazon.awssdk.services.s3.S3Client;
import software.amazon.awssdk.services.s3.S3ClientBuilder;
import software.amazon.awssdk.services.s3.model.GetObjectRequest;
import software.amazon.awssdk.services.s3.model.GetObjectResponse;
import software.amazon.awssdk.services.s3.model.PutObjectRequest;
import software.amazon.awssdk.services.s3.model.PutObjectResponse;

public class S3Presigner {

	private Region region;
	private AwsCredentialsProvider awsCredentialsProvider;
	private Duration expirationTime;
	private Integer timeOffset;
	private S3PresignExecutionInterceptor presignInterceptor;
	private S3Client s3Client;

	private S3Presigner() {
		presignInterceptor = new S3PresignExecutionInterceptor();
	}

	public static Builder builder() {
		return new Builder();
	}

	public URI presignS3DownloadLink(String bucketName, String fileName) throws SdkClientException {
		try {

			GetObjectRequest s3GetRequest = GetObjectRequest.builder().bucket(bucketName).key(fileName).build();
			ResponseInputStream<GetObjectResponse> response = s3Client.getObject(s3GetRequest);
			response.close();

			return presignInterceptor.getSignedURI();
		} catch (Throwable t) {
			if (t instanceof SdkClientException) {
				throw (SdkClientException) t;
			}
			throw SdkClientException.builder().cause(t).build();
		}
	}

	public URI presignS3UploadLink(String bucketName, String fileName) throws SdkClientException {
		try {

			PutObjectRequest s3PutRequest = PutObjectRequest.builder().bucket(bucketName).key(fileName).build();
			PutObjectResponse response = s3Client.putObject(s3PutRequest, RequestBody.empty());

			return presignInterceptor.getSignedURI();
		} catch (Throwable t) {
			if (t instanceof SdkClientException) {
				throw (SdkClientException) t;
			}
			throw SdkClientException.builder().cause(t).build();
		}
	}

	public static class Builder {
		S3Presigner presigner = new S3Presigner();

		public S3Presigner build() {
			if (presigner.awsCredentialsProvider == null) {
				DefaultCredentialsProvider provider = DefaultCredentialsProvider.create();
				presigner.awsCredentialsProvider = provider;
			}

			if (presigner.region == null) {
				presigner.region = new DefaultAwsRegionProviderChain().getRegion();
			}

			if (presigner.expirationTime == null) {
				presigner.expirationTime = Duration.ofDays(4);
			}

			if (presigner.timeOffset == null) {
				presigner.timeOffset = 2;
			}

			S3ClientBuilder s3Builder = S3Client.builder().region(presigner.region).credentialsProvider(presigner.awsCredentialsProvider);
			s3Builder.overrideConfiguration(ClientOverrideConfiguration.builder().addExecutionInterceptor(presigner.presignInterceptor).build());
			s3Builder.httpClient(new NullSdkHttpClient());
			presigner.s3Client = s3Builder.build();

			return presigner;
		}

		public Builder awsCredentials(AwsCredentialsProvider awsCredentialsProvider) {
			presigner.awsCredentialsProvider = awsCredentialsProvider;
			return this;
		}

		public Builder region(Region region) {
			presigner.region = region;
			return this;
		}

		public Builder expiration(Duration expirationTime) {
			presigner.expirationTime = expirationTime;
			return this;
		}

		public Builder timeOffset(Integer timeOffset) {
			presigner.timeOffset = timeOffset;
			return this;
		}

	}

	public static class NullSdkHttpClient implements SdkHttpClient {

		@Override
		public void close() {

		}

		@Override
		public ExecutableHttpRequest prepareRequest(HttpExecuteRequest request) {
			return new ExecutableHttpRequest() {
				@Override
				public HttpExecuteResponse call() {
					return HttpExecuteResponse.builder().response(SdkHttpResponse.builder().statusCode(200).build()).responseBody(AbortableInputStream.create(new ByteArrayInputStream(new byte[0]))).build();
				}

				@Override
				public void abort() {
				}
			};
		}
	}

	public class S3PresignExecutionInterceptor implements ExecutionInterceptor {

		final private AwsS3V4Signer signer;
		private URI signedURI;

		public S3PresignExecutionInterceptor() {
			this.signer = AwsS3V4Signer.create();
		}

		@Override
		public void beforeTransmission(BeforeTransmission context, ExecutionAttributes executionAttributes) {
			// remove all headers because a Browser that downloads the shared URL will not send the exact values. X-Amz-SignedHeaders should only contain the host header.
			SdkHttpFullRequest modifiedSdkRequest = (SdkHttpFullRequest) context.httpRequest().toBuilder().clearHeaders().build();


			executionAttributes.putAttribute(AwsSignerExecutionAttribute.AWS_CREDENTIALS, awsCredentialsProvider.resolveCredentials());
			executionAttributes.putAttribute(AwsSignerExecutionAttribute.PRESIGNER_EXPIRATION, Instant.ofEpochSecond(System.currentTimeMillis()/1000).plus(expirationTime));
			executionAttributes.putAttribute(AwsSignerExecutionAttribute.SERVICE_SIGNING_NAME, "s3");
			executionAttributes.putAttribute(AwsSignerExecutionAttribute.SIGNING_REGION, region);
			executionAttributes.putAttribute(AwsSignerExecutionAttribute.TIME_OFFSET, timeOffset);
			SdkHttpFullRequest signedRequest = signer.presign(modifiedSdkRequest, executionAttributes);// sign(getRequest, new ExecutionAttributes());
			signedURI = signedRequest.getUri();
		}

		public URI getSignedURI() {
			return signedURI;
		}

	}

}
