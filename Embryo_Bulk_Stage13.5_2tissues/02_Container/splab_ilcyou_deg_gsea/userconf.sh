#!/usr/bin/with-contenv bash

## Set defaults for environmental variables in case they are undefined
USER=${USER:=rstudio}
PASSWORD=${PASSWORD:=rstudio}
USERID=${USERID:=10000}
GROUPID=${GROUPID:=10000}
ROOT=${ROOT:=FALSE}
UMASK=${UMASK:=022}

touch /tmp/log_user.txt

echo "USER=$USER" >> /tmp/log_user.txt
echo "PASSWORD=$PASSWORD" >> /tmp/log_user.txt
echo "USERID=$USERID" >> /tmp/log_user.txt
echo "GROUPID=$GROUPID" >> /tmp/log_user.txt
echo "ROOT=$ROOT" >> /tmp/log_user.txt
echo "UMASK=$UMASK" >> /tmp/log_user.txt

if [ "$USERID" -lt 1000 ]
# Probably a macOS user, https://github.com/rocker-org/rocker/issues/205
  then
    echo "$USERID is less than 1000, setting minumum authorised user to 499" >> /tmp/log_user.txt
    echo auth-minimum-user-id=499 >> /etc/rstudio/rserver.conf >> /tmp/log_user.txt
fi

if [ "$USERID" -ne 999 ]
## Configure user with a different USERID if requested.
  then
    echo "deleting user rstudio" >> /tmp/log_user.txt
    userdel rstudio
    echo "creating new $USER with UID $USERID" >> /tmp/log_user.txt
    useradd -m $USER -u $USERID
    mkdir /home/$USER
    chown -R $USER /home/$USER
    usermod -a -G staff $USER
elif [ "$USER" != "rstudio" ]
  then
    ## cannot move home folder when it's a shared volume, have to copy and change permissions instead
    cp -r /home/rstudio /home/$USER
    ## RENAME the user   
    usermod -l $USER -d /home/$USER rstudio
    groupmod -n $USER rstudio
    usermod -a -G staff $USER
    chown -R $USER:$USER /home/$USER
    echo "USER is now $USER"   >> /tmp/log_user.txt
fi

if [ "$GROUPID" -ne 999 ]
## Configure the primary GID (whether rstudio or $USER) with a different GROUPID if requested.
  then
    echo "Modifying primary group $(id $USER -g -n)"  >> /tmp/log_user.txt
    groupmod -g $GROUPID $(id $USER -g -n)
    echo "Primary group ID is now custom_group $GROUPID"  >> /tmp/log_user.txt
fi
  
## Add a password to user
echo "$USER:$PASSWORD" | chpasswd

# Use Env flag to know if user should be added to sudoers
if [ "$ROOT" == "TRUE" ]
  then
    adduser $USER sudo && echo '%sudo ALL=(ALL) NOPASSWD:ALL' 
    echo "$USER added to sudoers"  >> /tmp/log_user.txt
fi

## Change Umask value if desired
if [ "$UMASK" -ne 022 ]
  then
    echo "server-set-umask=false" >> /etc/rstudio/rserver.conf
    echo "Sys.umask(mode=$UMASK)" >> /home/$USER/.Rprofile
fi

## add these to the global environment so they are avialable to the RStudio user 
echo "HTTR_LOCALHOST=$HTTR_LOCALHOST" >> /etc/R/Renviron.site
echo "HTTR_PORT=$HTTR_PORT" >> /etc/R/Renviron.site
