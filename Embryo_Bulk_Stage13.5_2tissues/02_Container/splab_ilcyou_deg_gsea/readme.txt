
This image contains:

 - R 3.4.1
 - Rstudio server (installation requires the userconf.sh file)
 - Deseq2
  

# ######################
     COMPILE THE IMAGE
# ######################

docker build -t splab_ilcyou_deg_gsea ~/workspace/ciml-docker/images/project/SPlab/ILCyou/deg_gsea

# ######################
     RUN THE IMAGE
# ######################

docker run -d --name splab_ilcyou_deg_gsea -p 9494:8787 -v /mnt:/mnt -v /home/$USER/workspace:/home/$USER/workspace -e PASSWORD= -e USER=$(whoami) -e USERID=$(id -u) -e GROUPID=$(id -g) splab_ilcyou_deg_gsea
 
# ######################
     CONNECT TO RSTUDIO
# ######################
 
 In an Internet browser, type as url : http://127.0.0.1:8787 and use the login/password: <USER>/rstudio
 
# ######################
	 NOTES
# ######################
 
 - To use knitr PDF compilation instead of Sweave, you have to go into Rstudio menu Tools->Global Options->Sweave->Weave Rnw files with.. and select "knitr".
 
 
 