###
# Do some basic differential expression analysis on the expression data.
#

require 'catpaws'

#generic settings
set :aws_access_key,  ENV['AMAZON_ACCESS_KEY']
set :aws_secret_access_key , ENV['AMAZON_SECRET_ACCESS_KEY']
set :ec2_url, ENV['EC2_URL']
set :ssh_options, { :user => "ubuntu", :keys=>[ENV['EC2_KEYFILE']]}
set :key, ENV['EC2_KEY']
set :key_file, ENV['EC2_KEYFILE']
set :ami, 'ami-52794c26'  #EC2 eu-west-1 32bit Lucid
set :instance_type, 'm1.small'
set :s3cfg, ENV['S3CFG'] #location of ubuntu s3cfg file
set :working_dir, '/mnt/work'

set :nhosts, 1
set :group_name, 'ctx_bmp_lumiarray'

set :snap_id, `cat SNAPID`.chomp #ec2 eu-west-1 
set :vol_id, `cat VOLUMEID`.chomp #empty until you've created a new volume
set :ebs_size, 5  #give us enough space to generate some analysis data
set :availability_zone, 'eu-west-1a'  #wherever the ami is
set :dev, '/dev/sdf'
set :mount_point, '/mnt/data'



#make a new EBS volume from this snap 
#cap EBS:create

#and mount your EBS
#cap EBS:attach
#cap EBS:format_xfs
#cap EBS:mount_xfs



desc "install R on all running instances in group group_name"
task :install_r, :roles  => group_name do
  user = variables[:ssh_options][:user]
  sudo 'apt-get update'
  sudo 'apt-get -y install r-base'
  sudo 'apt-get -y install build-essential libxml2 libxml2-dev libcurl3 libcurl4-openssl-dev'
  upload("scripts/R_setup.R", "#{working_dir}/R_setup.R")
  run "cd #{working_dir} && chmod +x R_setup.R"
  sudo "Rscript #{working_dir}/R_setup.R"
end
before "install_r", "EC2:start"
  


desc "run pre-processing on expression data"
task :pp_expression_data, :roles => group_name do
  run "mkdir -p #{working_dir}/scripts"
  upload "scripts/limma_xpn.R", "#{working_dir}/scripts/limma_xpn.R"
  run "chmod +x #{working_dir}/scripts/limma_xpn.R"
  run "cd #{mount_point} && Rscript #{working_dir}/scripts/limma_xpn.R #{mount_point}/Non-background_subtracted_Sample_Probe_Profile.txt #{mount_point}/limma_results.csv"
end
before "pp_expression_data", "EC2:start"


desc "run QC checks on the pre-processed quality control"
task "pp_qc_expression_data", :roles => group_name do
  puts "TODO"
end  
before "pp_qc_exprssion_data", "EC2:start"


desc "Fetch ReMoat data which has mm9 probe positions"
task :get_remoat_anno, :roles => group_name do
  sudo "apt-get -y install unzip"
  run "mkdir -p #{working_dir}/lib"
  run "rm -Rf  #{working_dir}/lib/Annotation_Illumina_Mouse*"
  run "cd #{working_dir}/lib && curl http://www.compbio.group.cam.ac.uk/Resources/Annotation/final/Annotation_Illumina_Mouse-WG-V1_mm9_V1.0.0_Aug09.zip > Annotation_Illumina_Mouse-WG-V1_mm9_V1.0.0_Aug09.zip "
  run "cd #{working_dir}/lib && unzip Annotation_Illumina_Mouse-WG-V1_mm9_V1.0.0_Aug09.zip"
end
before 'get_remoat_anno', 'EC2:start'



desc "Make an IRanges RangedData object from expression data"
task :xpn2rd, :roles => group_name do
  user = variables[:ssh_options][:user]
  run "mkdir -p #{working_dir}/scripts"
  upload("scripts/xpn_csv_to_iranges.R", "#{working_dir}/scripts/xpn_csv_to_iranges.R")
  run "cd #{working_dir}/scripts && chmod +x xpn_csv_to_iranges.R"
  run "Rscript #{working_dir}/scripts/xpn_csv_to_iranges.R #{mount_point}/limma_results.csv #{working_dir}/lib/Annotation_Illumina_Mouse-WG-V1_mm9_V1.0.0_Aug09.txt"

end
before "xpn2rd","EC2:start"


desc "Fetch expression data results"
task :get_xpn, :roles=> group_name do
  `mkdir -p results`
  download("#{mount_point}/limma_rd.csv", "results/limma_rd.csv")
end
before "get_xpn", "EC2:start"


#if you want to keep the results

#cap EBS:snapshot

#and then shut everything down:

# cap EBS:unmount
# cap EBS:detach
# cap EBS:delete - unless you're planning to use it again.
# cap EC2:stop




