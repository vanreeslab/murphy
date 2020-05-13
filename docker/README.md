### Gitlab
```
docker build -f Dockerfile_gitlab -t immc/murphy-ci:v1.1 .
docker login
docker push immc/murphy-ci:v1.0
```
to test the image locally:
```
docker run -it --rm immc/murphy-ci:v1.0
```


### GitHub
```
docker build -t immc/murphy:v1.1 .
docker login
docker push immc/murphy:v1.1
```
