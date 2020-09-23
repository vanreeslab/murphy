### GitLab
```
# download flups locally
git clone git@git.immc.ucl.ac.be:examples/flups.git 
# build the container
docker build -t immc/murphy-ci:v1.2 .
# push it
docker login
docker push immc/murphy-ci:v1.2
rm -rf flups
```
to test the image locally:
```
docker run -it --rm immc/murphy-ci:v1.2
```


### GitHub
```
#download flups locally
git clone git@git.immc.ucl.ac.be:examples/flups.git
# build the container
docker build --no-cache -t vanreeslab/murphy:v1.5 .
#push it
docker login
docker push vanreeslab/murphy:v1.5
```
to test the container locally:
```
docker run --rm -it vanreeslab/murphy:v1.5
```
