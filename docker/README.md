# build image
```
docker build -t immc/murphy-ci:v1.0 .
```

# push image to docker-hub
```
docker login
docker push immc/murphy-ci:v1.0
```

# local test:
```
docker run -it --rm immc/murphy-ci:v1.0
```
