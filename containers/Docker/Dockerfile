FROM sarusso/minimalmetadesktop:v0.2.0
MAINTAINER Stefano Alberto Russo <stefano.russo@gmail.com>

# Switch to root user
USER root

# Install APT dependencies
RUN apt-get install -y python3 python3-pip python3-wxgtk4.0 python3-wxgtk-webview4.0 python3-wxgtk-media4.0

# Install PIP dependencies
RUN pip3 install astropy==4.0.0 \
                 scipy==1.3.1 \
                 specutils==0.6 \
                 numpy==1.17.3 \
                 lmfit==0.9.14 \
                 cycler==0.10.0 \
                 statsmodels==0.10.1 \
                 matplotlib==3.1.1 \
                 sphinx==2.2.0\
                 tqdm==4.40.1

# Astrocook X environment startup
COPY containers/Docker/files/xstartup /opt/tigervnc/
RUN chmod 755 /opt/tigervnc/xstartup

# Prepare file for logs
RUN touch /var/log/astrocook.log && chown metauser:metauser  /var/log/astrocook.log

# Switch back to metauser user
USER metauser

# Add Astrocook code (as metauser)
COPY ./ /metauser/astrocook
