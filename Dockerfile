# Use an official Python runtime as a parent image
FROM homebrew/brew:latest

# Set the working directory in the container to /app
WORKDIR /home/linuxbrew/topospam

# Add the current directory contents into the container at /app
COPY --chown=linuxbrew:linuxbrew . /home/linuxbrew/topospam

USER root
# Install bash if not already installed
RUN apt-get update && apt-get install -y bash python3-venv libx11-dev libgl-dev libxrender-dev python3-pip

USER linuxbrew
# Create the virtual environment directly in the target directory
RUN python3 -m venv /home/linuxbrew/topospam/topospam_env

# Activate the virtual environment and install dependencies
RUN /bin/bash -c "source /home/linuxbrew/topospam/topospam_env/bin/activate && pip install ."