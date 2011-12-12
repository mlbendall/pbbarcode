all: build install

build:
	python setup.py build

install:
	python setup.py install

develop:
	python setup.py develop

clean:
	rm -rf build/;\
	find . -name "*.egg-info" | xargs rm -rf;\
	find . -name "*.pyc" | xargs rm -rf;\
	rm -rf dist/

test:
	nosetests -v tests/* 
