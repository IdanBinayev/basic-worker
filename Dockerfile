# Use an official Python runtime as a parent image
FROM python:latest

# Set the working directory in the container
WORKDIR /tmp

# Copy the current directory contents into the container at /usr/src/app
COPY xgb_model.pkl .
COPY . .


# Install any needed packages specified in requirements.txt

#RUN locale-gen en_US.UTF-8
#RUN update-locale en_US.UTF-8
#ENV LC_ALL="en_US.utf8"

COPY requirements*.txt /tmp/
# install pip and python

RUN apt-get update
# RUN add-apt-repository ppa:deadsnakes/ppa -y

RUN pip3 install --upgrade setuptools
RUN pip3 install --upgrade pip
RUN pip3 install --upgrade pillow
RUN apt install python3-gi python3-gi-cairo gir1.2-gtk-3.0 -y
#RUN pip3 install pygobject
RUN pip install --upgrade cython
RUN pip install wheel
RUN pip install Cmake
RUN pip3 install -r requirements1.txt
# RUN pip3 install -r requirements2.txt
# RUN pip3 install -r requirements3.txt
# RUN pip3 install -r requirements4.txt
# RUN pip3 install -r requirements5.txt
# RUN pip3 install -r requirements6.txt
# RUN pip3 install -r requirements7.txt
# RUN pip3 install -r requirements8.txt --ignore-installed
RUN pip3 install boto3
RUN apt-get update && apt-get install ffmpeg libsm6 libxext6  -y
RUN apt autoremove -y
COPY step_py_runner /srv/scripter/step_py_runner
WORKDIR /srv/scripter
COPY . /srv/scripter
