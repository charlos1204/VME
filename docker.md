# Docker installation for running the container (optional)
Only for running the Docker image of the code.
* to install Docker in Debian based distribution (Ubuntu):<br>
`sudo apt update`<br>
`sudo apt-get install apt-transport-https ca-certificates curl gnupg-agent software-properties-common`<br>
`curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -`<br>
`sudo apt-key fingerprint 0EBFCD88`<br>
`sudo add-apt-repository "deb [arch=amd64] https://download.docker.com/linux/ubuntu $(lsb_release -cs) stable"`<br>
`sudo apt-get update`<br>
`sudo apt-get install docker-ce docker-ce-cli containerd.io`

* To install Docker in Windows 10 follow the steps [here](https://docs.docker.com/docker-for-windows/install/)
* To install Docker in macOS follow the steps [here](https://docs.docker.com/docker-for-mac/install/)

Run this example in a terminal to make sure Docker is installed correctly:<br>
`docker run hello-world`

**NOTE**<br>
In linux you have to add sudo before docker!
