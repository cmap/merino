import boto3
import botocore

def stream_file_from_s3(bucket_name, file_key):
    s3 = boto3.resource("s3")
    file =s3.Object(bucket_name, file_key)
    file_content = file.get()['Body'].read()
    return file_content

def download_file_from_s3(bucket_name, file_key, local_path):
    s3=boto3.resource('s3')
    try:
        print 'Reading file {} from bucket {}'.format(file_key, bucket_name)
        s3.Bucket(bucket_name).download_file(file_key, local_path)

    except botocore.exceptions.ClientError as e:

        if e.response['Error']['Code'] == "404":
            error_message = "The file located at {} from bucket {} does not exist".format(file_key, bucket_name)
            print error_message

        else:
            error_message = "Failed to download file located at {} from bucket {}".format(file_key, bucket_name)
            print error_message

        raise Exception(e)

def upload_file_to_s3(local_file, bucket_name, file_key):
    s3 = boto3.client('s3')
    try:
        s3.upload_fileobj(local_file, Bucket=bucket_name, Key=file_key)
    except boto3.exceptions.S3UploadFailedError as error:
        print "file {} failed to upload to s3 bucket {} at key {} because {}".format(local_file, bucket_name, file_key, error)
        raise Exception(error)