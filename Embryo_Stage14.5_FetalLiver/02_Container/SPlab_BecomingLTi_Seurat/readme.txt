This image contains:

 - R 3.6.1
 - Rstudio server
 - Packages for the Seurat analysis
  

# ######################
     COMPILE THE IMAGE
# ######################

docker build -t splab_becominglti_seurat /mnt/NAS7/Workspace/spinellil/ciml-bip/project/SPlab/BecomingLTi/Embryo_Stage13.5_FetalLiver/00_Container/SPlab_BecomingLTi_Seurat

# ######################
     RUN THE IMAGE
# ######################

docker run --name splab_becominglti_seurat -d -p 8787:8787 -v /mnt:/mnt -e PASSWORD= -e USER=$(whoami) -e USERID=$(id -u) -e GROUPID=$(id -g)  splab_becominglti_seurat
 
# ######################
     CONNECT TO RSTUDIO
# ######################
 
 In an Internet browser, type as url : http://127.0.0.1:8787 and use the login/password: <user>/rstudio
 
 
