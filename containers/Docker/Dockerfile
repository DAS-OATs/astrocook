FROM ubuntu:18.04
MAINTAINER Stefano Alberto Russo <stefano.russo@gmail.com>

#----------------------
# Basics
#----------------------

# Set non-interactive
ENV DEBIAN_FRONTEND noninteractive

# Update first of all
RUN apt-get update

# Utilities
RUN apt-get install -y nano telnet unzip wget supervisor build-essential python-dev git-core openjdk-8-jre


#------------------------
# "Meta" user
#------------------------

# Add group. We chose GID 65527 to try avoiding conflicts.
RUN groupadd -g 65527 metauser

# Add user. We chose UID 65527 to try avoiding conflicts.
RUN useradd metauser -d /metauser -u 65527 -g 65527 -m -s /bin/bash

# Add metuaser user to sudoers
RUN adduser metauser sudo

# Install suodo
RUN apt-get install sudo -y

# No pass sudo (for everyone, actually)
COPY containers/Docker/files/sudoers /etc/sudoers

#------------------------
# Supervisor conf
#------------------------

# Supervisord conf
COPY containers/Docker/files/supervisord.conf /etc/supervisor/

#------------------------
# VNC
#------------------------

# Install xvfb that triggers minimal install of X base packages and xterm as sample application
RUN apt-get install xvfb xterm  -y

# Install base packages for VNC server and headless desktop (2)
RUN cd /opt && wget https://bintray.com/tigervnc/stable/download_file?file_path=tigervnc-1.8.0.x86_64.tar.gz -O tigervnc-1.8.0.x86_64.tar.gz \
            && tar -zxvf tigervnc-1.8.0.x86_64.tar.gz \
            && mv tigervnc-1.8.0.x86_64 tigervnc

# Supervisord configuration
COPY containers/Docker/files/supervisord_vnc.conf /etc/supervisor/conf.d/
COPY containers/Docker/files/run_vnc.sh /etc/supervisor/conf.d/
RUN chmod 755 /etc/supervisor/conf.d/run_vnc.sh

# X startup
COPY containers/Docker/files/basic_xstartup /opt/tigervnc/xstartup
RUN chmod 755 /opt/tigervnc/xstartup


# Web VNC (noVNC) v0.6.1.
# NOTE: this is a custom version from Doro Wu (fcwu.tw@gmail.com).
# TODO: Check differences and maybe move to 0.6.2
COPY containers/Docker/files/noVNC.tar.gz /usr/lib/
RUN cd /usr/lib/ && tar -zxvf noVNC.tar.gz
COPY containers/Docker/files/index.html /usr/lib/noVNC

RUN apt-get install -y net-tools

# Supervisord configuration
COPY containers/Docker/files/supervisord_webvnc.conf /etc/supervisor/conf.d/

# Install APT dependencies
RUN apt-get install -y python3 python3-pip python3-wxgtk4.0 python3-wxgtk-webview4.0 python3-wxgtk-media4.0

# Install PIP dependencies
RUN pip3 install astropy==3.2.2 \
                 scipy==1.3.1 \
                 specutils==0.6 \
                 numpy==1.17.3 \
                 lmfit==0.9.14 \
                 cycler==0.10.0 \
                 statsmodels==0.10.1 \
                 matplotlib==3.1.1 \
                 sphinx==2.2.0

# Astroccok X environment setup/startup
RUN apt-get install fluxbox -y
COPY containers/Docker/files/xstartup /opt/tigervnc/
RUN chmod 755 /opt/tigervnc/xstartup

COPY containers/Docker/files/dot_fluxbox /metauser/.fluxbox
RUN chown -R metauser:metauser /metauser/.fluxbox
COPY containers/Docker/files/cookies.jpg /usr/share/images/fluxbox/cookies.jpg

# Prepare file for logs
RUN touch /var/log/astrocook.log && chown metauser:metauser  /var/log/astrocook.log

# Add Astrocook code (as metauser) adn switch back to root
USER metauser
COPY ./ /metauser/astrocook
USER root

#----------------------
# Entrypoint
#----------------------

# Copy entrypoint
COPY containers/Docker/files/entrypoint.sh /

# Give right permissions
RUN chmod 755 /entrypoint.sh

# Set entrypoint
ENTRYPOINT ["/entrypoint.sh"]

# To access: expose 8590/tcp and 5900/tcp
