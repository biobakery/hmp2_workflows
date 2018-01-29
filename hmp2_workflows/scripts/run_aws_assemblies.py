# -*- coding: utf-8 -*-

"""
run_aws_assemblies.py
~~~~~~~~~~~~~~~~~~~~~

Setups and fans out assemblies of the HMP2 MGX dataset on AWS.

Copyright (c) 2017 Harvard School of Public Health

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
    THE SOFTWARE.
"""


import argparse
import os
import sys
import tempfile

from time import sleep

import boto3

from botocore.exceptions import ClientError
from paramiko import SSHClient, AutoAddPolicy
from scp import SCPClient


def parse_cli_arguments():
    """
    Parses any command-line arguments passed into this script.

    Args:
        None

    Requires:
        None

    Returns:
        argparse.ArgumentParser: Object containing all command-line args.
    """
    parser = argparse.ArgumentParser('Sets up and fans out HMP2 metagenomic '
                                     'assemblies on AWS.')
    parser.add_argument('-f', '--forward-read', required=True,
                        help='Forward metagenomic sequence read.')
    parser.add_argument('-r', '--reverse-read', required=True,
                        help='Reverse metagenomic sequence read.')
    parser.add_argument('-o', '--output-dir', required=True,
                        help='Output directory to download assembled '
                        'metagenome too.')
    parser.add_argument('-aws_a', '--access-key', required=True,
                        help='AWS access key.')
    parser.add_argument('-aws_s', '--secret-key', required=True,
                        help='AWS secret key.')                        
    parser.add_argument('-k', '--ssh-key', required=True,
                        help='SSH key used to connect to EC2 instances.') 
    parser.add_argument('-ami', '--ami-id', default='ami-14bed66e',
                        help='OPTIONAL. AMI ID for AWS image containing '
                        'IGS assembly pipeline.')

    return parser.parse_args()    


def generate_mapping_file(f_read, r_read, out_file):
    """Generates the required IGS assembly pipeline mapping file.

    Args:
        f_read (string): Path to forward sequence read.
        r_read (string): Path to reverse sequence read.
        out_file (string): Path to output mapping file.

    Requires:
        None

    Returns: 
        Path to output mapping file
    """
    f_basename = os.path.basename(f_read)
    r_basename = os.path.basename(r_read)
    sample_base = f_basename.split(os.extsep)[0].replace('_R1', '')

    out_file.write("%s\t1\t%s\n" % (sample_base, 
                                    os.path.join('/mnt', 'data', f_basename.replace('.gz', ''))))
    out_file.write("%s\t2\t%s\n" % (sample_base,
                                    os.path.join('/mnt', 'data', r_basename.replace('.gz', ''))))
    out_file.close()

    return out_file.name


def upload_and_process_input_files(ssh, instance, ssh_key, input_files):
    """Process input files to generate a mapping file for the IGS assembly 
    pipeline and uploads the provided sequence files (paired-end) along with 
    the mapping file to the instance.

    Args:
        ssh (paramiko.SSHClient): Paramiko SSH client used to interface
            with the EC2 instance.
        instance (boto3.EC2.Instance): boto3 Instance representation of EC2
            instance.
        ssh_key (string): Path to the public SSH keypair used to connect to the
            EC2 instance.
        input_files (list): The paired-end sequences to upload to the instance.

    Requires:
        None

    Returns:
        paramiko.SSHClient: The paramiko SSH connection to our EC2 instance.
        scp.SCPClient: An SCP client connection to our EC2 instance.
    """
    ssh.connect(hostname=instance.public_dns_name, 
                username='ec2-user',
                key_filename=ssh_key)

    mapping_file = generate_mapping_file(input_files[0],
                                         input_files[1],
                                         tempfile.NamedTemporaryFile(delete=False))

    scp = SCPClient(ssh.get_transport())

    print "Transferring file %s to EC2 instance..." % input_files[0],
    sys.stdout.flush()
    scp.put(input_files[0], remote_path='/mnt/data/')
    print "    DONE"

    print "Transferring file %s to EC2 instance..." % input_files[1],
    sys.stdout.flush()
    scp.put(input_files[1], remote_path='/mnt/data/')
    print "    DONE"

    print "Transferring file %s to EC2 instance..." % mapping_file,
    sys.stdout.flush()
    scp.put(mapping_file, remote_path='/mnt/data/zz00_input_locations.txt')
    print "    DONE"

    print "Extracting compressed sequence files...",
    sys.stdout.flush()
    (stdin, stdout, stderr) = ssh.exec_command('gunzip /mnt/data/*.gz')
    stdout.channel.recv_exit_status() ## Block until done
    print "    DONE"

    print "Running IMA_setup...",
    sys.stdout.flush()
    (stdin, stdout, stderr) = ssh.exec_command('cd /mnt/data; '
                                               '/bin/bash -c /home/ec2-user/bin/IMA_setup')
    stdout.channel.recv_exit_status() ## Block again!
    print "    DONE"

    os.remove(mapping_file)

    return scp


def start_ec2_instance(ec2r, ec2c, ami_id):
    """Starts up an AWS instance containing the IGS assembly pipeline
    and configures the instance so that it is ready to run on HMP2 
    metagenomics data.

    Args:
        ec2r (boto3.resource): The boto3 AWS EC2 resource interface.
        ec2c (boto3.client): The boto3 AWS EC2 client interface.
        ami_id (string): The AMI ID for the AWS image containing the requisite
            software.

    Requires:
        None

    Returns:
        The instance ID for the instance started                 

    """
    userdata = """#cloud-config
        runcmd:
            - [ sh, -c, "parted -s -a optimal /dev/nvme0n1 mklabel msdos mkpart primary 0% 100%" ]
            - [ sh, -c, "mkfs /dev/nvme0n1p1" ]
            - [ sh, -c, "mkdir -p /mnt/data" ]
            - [ sh, -c, "mount /dev/nvme0n1p1 /mnt/data" ]
            - [ sh, -c, "chmod ugo+rwx /mnt/data" ]
    """ 
    waiter = ec2c.get_waiter('instance_running')
 
    try:
        print "Starting EC2 instance...",
        sys.stdout.flush()
        response = ec2r.create_instances(ImageId=ami_id,
                                         InstanceType='i3.2xlarge',
                                         KeyName='hmp2_keypair',
                                         SecurityGroupIds=['sg-dc6f67a9'],
                                         EbsOptimized=True,
                                         UserData=userdata,
                                         MinCount=1,
                                         MaxCount=1,
                                         DryRun=False)
        instance_id = response[0].id
        waiter.wait(InstanceIds=[instance_id])
        sleep(100)
        print "    DONE"
    except ClientError as e:
        print e        

    instance = ec2r.Instance(instance_id)
    return instance


def start_assembly_pipeline(ssh):
    """Starts up the IGS assembly pipeline on the provided EC2 instance.

    Args:
        ssh (paramiko.SSHClient): Paramiko SSH connection to the assembly EC2
            instance

    Requires:
        None

    Returns:
        int: The process ID for the pipeline run   
    """
    print "Starting assembly pipeline...",
    sys.stdout.flush()
    (stdin, stdout, stderr) = ssh.exec_command('cd /mnt/data; '
                                               'nohup /bin/bash '
                                               '-lc \'/mnt/data/003_run_assemblies.sh '
                                               '> /dev/null 2>&1 &\'')
    stdout.channel.recv_exit_status() ## Block again!
    print "    DONE"

    (stdin, stdout, stderr) = ssh.exec_command('pgrep -f 003_run')
    pid = int(stdout.readline())
 

    return pid


def is_assembly_running(ssh):
    """Checks to see if IGS assembly pipeline is still running on the provided
    EC2 instance.

    Args:
        ssh (paramiko.SSHClient): Paramiko SSH connection to the assembly EC2 
            instance.

    Requires:
        None

    Returns:
        boolean: True if assembly pipeline is still running otherwise False.                
    """
    is_running = True
    (stdin, stdout, stderr) = ssh.exec_command('pgrep -f 003_run')
    exit_status = stdout.channel.recv_exit_status()

    return True if exit_status == 0 else False


def download_assembly_files(ssh, scp, sample_base, output_dir):
    """Downloads the complete assembly files from the EC2 instance.

    Args:
        ssh (paramiko.SSHClient): The paramiko SSH connection to the assembly
            instance.
        scp (scp.SCPClient): An SCP client/connection to the assembly instance.
        sample_base (string): The sample name for which assembly is being run
            on.
        output_dir (string): The path to the desired location to download all 
            assembly files too.

    Requires:
        None

    Returns:
        list: All files downloaded from the EC2 instance.                    

    """
    scp.get('/mnt/data/%s/%s__FINAL_ASSEMBLY_UNIQUE_IDS.consolidated.fna' % (sample_base, sample_base), 
            local_path=output_dir)
    scp.get('/mnt/data/002_assembly_logs/',
            recursive=True,
            local_path=output_dir)
    scp.get('/mnt/data/%s/01_logs/' % sample_base, 
            recursive=True,
            local_path=output_dir)            
    scp.get('/mnt/data/%s/07_consolidated_assembly/999_FINAL_OUTPUT/' % sample_base,
            recursive=True,
            local_path=output_dir)


def main(args):
    input_files = [args.forward_read, args.reverse_read]
    sample_base = (os.path.basename(input_files[0])
                     .split(os.extsep)[0]
                     .replace('_R1', ''))
    output_dir = os.path.join(args.output_dir, sample_base)

    try:
       os.mkdir(output_dir)
    except OSError, e:
        if e.errno != os.errno.EEXIST:
            raise
        pass       

    ec2r = boto3.resource('ec2',
                          aws_access_key_id=args.access_key,
                          aws_secret_access_key=args.secret_key,
                          region_name='us-east-1')
    ec2c = boto3.client('ec2',
                         aws_access_key_id=args.access_key,
                         aws_secret_access_key=args.secret_key,
                         region_name='us-east-1')
    ssh_client = SSHClient()
    ssh_client.set_missing_host_key_policy(AutoAddPolicy())

    try:
        instance = start_ec2_instance(ec2r, ec2c, args.ami_id)
        scp_client = upload_and_process_input_files(ssh_client, instance, args.ssh_key, input_files)
        assembly_pid = start_assembly_pipeline(ssh_client)

        print "Running assembly...",
        sys.stdout.flush()
        is_running = True
        while is_running:
            sleep(180)
            is_running = is_assembly_running(ssh_client)

        print "    DONE"

        print "Downloading assembly files...",
        sys.stdout.flush()
        download_assembly_files(ssh_client, scp_client, sample_base, output_dir)
        print "    DONE"

        ## Shut-down the instance once we are done here
        instance.terminate()
    except Exception:
        instance.terminate()        
        raise
    finally:
        instance.terminate()        


if __name__ == "__main__":
    main(parse_cli_arguments())