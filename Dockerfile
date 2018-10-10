FROM ubuntu:16.04

MAINTAINER	Andrey Smelter

# Copy ssc code to the root directory and set working directory
COPY . /ssc/
WORKDIR /ssc

# Install gcc and python
RUN apt-get update
RUN apt-get install gcc -y
RUN apt-get install python3 -y
RUN apt-get install python3-pip -y

# Install ssc python requirements
RUN pip3 install -r /ssc/requirements.txt

# Change registration algorithm executable permissions
RUN chmod a+x /ssc/ssc/bin/CRS_EXE

# Set entry point
ENTRYPOINT [ "python3", "-m", "ssc" ]
