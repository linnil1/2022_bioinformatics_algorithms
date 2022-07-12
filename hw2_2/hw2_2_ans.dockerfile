from docker.io/library/python:3.10
copy HW2.2/Answer/pyramid.py /pyramid.py
entrypoint ["python3", "/pyramid.py"]
cmd ["10"]
