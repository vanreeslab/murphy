### GitLab
```
# download flups locally
git clone git@git.immc.ucl.ac.be:examples/flups.git 
# build the container
docker build -t immc/murphy-ci:v1.6 .
# push it
docker login
docker push immc/murphy-ci:v1.6
rm -rf flups
```
to test the image locally:
```
docker run -it --rm immc/murphy-ci:v1.6
```


### GitHub
```
#download flups locally
git clone git@git.immc.ucl.ac.be:examples/flups.git
# build the container
docker build --no-cache -t vanreeslab/murphy:v1.8 .
#push it
docker login
docker push vanreeslab/murphy:v1.8
```
to test the container locally:
```
docker run --rm -it vanreeslab/murphy:v1.8
```
or to permanently build it and mount the murphy folder
```
docker run -it --name murphy -v /Users/tgillis/Dropbox/research/codes/murphy:/murphy vanreeslab/murphy:v1.8
```
