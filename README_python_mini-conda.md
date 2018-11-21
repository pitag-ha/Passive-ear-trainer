# Miniconda on Heroku Example App

This repository contains two things:

- A `Dockerfile`, which installs [scikit-learn](http://scikit-learn.org/stable/) with [miniconda](http://conda.pydata.org/miniconda.html), and a few [pip](https://pip.pypa.io/en/stable/) dependencies.
- A [Flask](http://flask.pocoo.org) `webapp`, which utilizes basic functionality of `scikit-learn`.

All [Anaconda packages](https://docs.continuum.io/anaconda/pkg-docs) are supported—`scikit-learn` is just being used here as an example. 

## ☤ Advantages over [Conda Buildpack](https://github.com/kennethreitz/conda-buildpack):

- No slug size limit (Anaconda packages can be very large). 
- Exact Miniconda environment, from Continuum Analytics.

## ☤ Deploy this Application:

Deploy with the [Container Registry and Runtime](https://devcenter.heroku.com/articles/container-registry-and-runtime):
     $ sudo usermod -a -G docker $USER
     ($ heroku plugins:install heroku-container-registry) not working anymore...
     $ heroku plugins:install @heroku-cli/plugin-container-registry
     $ heroku container:login
     
     $ git clone https://github.com/heroku-examples/python-miniconda
     $ cd python-miniconda
     
     $ heroku create
     $ heroku container:push web
     ($ heroku:release web) also not working for me
     $ heroku container:release web
     $ heroku open 


When you're running out of disk:
Docker provides a single command that will clean up any resources — images, containers, volumes, and networks — that are dangling (not associated with a container):

    docker system prune

To additionally remove any stopped containers and all unused images (not just dangling images), add the -a flag to the command:

    docker system prune -a

docker images -a | grep "pattern" | awk '{print $3}' | xargs docker rmi

All the Docker images on a system can be listed by adding -a to the docker images command. Once you're sure you want to delete them all, you can add the -q flag to pass the Image ID to docker rmi:

List:

    docker images -a

Remove:

    docker rmi $(docker images -a -q)


✨🍰✨
