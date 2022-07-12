#/bin/bash
set -e -o pipefail

# Copy question folder
cp /home/weber50432/0712linux_class/HW2.2_package.tar.gz .
tar xvf HW2.2_package.tar.gz
cat HW2.2/Introduction.txt
mkdir -p answer

# Q1: find findme.txt in the HW_find_file directory
find HW2.2/HW_find_file/ -name "findme.txt" > answer/location.txt

# Q2: List owner_name and dir_name in /home
ls /home/ -alh | awk '{print $3 " " $9}' | grep -v "root" > answer/home.txt

# Q3: Build the docker
cp hw2_2_ans.dockerfile answer/Dockerfile
podman build . -f hw2_2_ans.dockerfile -t linnil1/summer2022_hw2_2
podman run -it --rm linnil1/summer2022_hw2_2

# zip all data
tar -czf hw2_2_ans.tar.gz answer
mv hw2_2_ans.tar.gz /home/weber50432/0712linux_class/Assignments/r08631020.tar.gz
