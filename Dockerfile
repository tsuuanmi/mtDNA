FROM ubuntu:22.04

# Install system dependencies
RUN apt update \
    && apt install -y build-essential python3 python3-pip parallel unzip \
    file libgl1-mesa-glx \
    gcc-multilib libc6-i386 lib32stdc++6

# Create a symbolic link so 'python' works as a command
RUN ln -sf /usr/bin/python3 /usr/bin/python

# Create app directory structure
WORKDIR /app

# Copy project files
COPY . /app

# Create directory structure as in setup.sh
RUN bash scripts/setup.sh

# Make the tools executable
RUN chmod -R +x /app/tools/*

# Install Python dependencies 
RUN pip3 install --no-cache-dir -r requirements.txt