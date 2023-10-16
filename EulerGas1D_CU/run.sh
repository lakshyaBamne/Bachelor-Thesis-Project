# make a result directory if it does not already exist
mkdir -p env
mkdir -p result
mkdir -p result1
mkdir -p result2
mkdir -p plots

# remove all the existing data files to create new ones
find . -name "*.txt" -type f -delete
find . -name "*.gif" -type f -delete
# find . -name "*.png" -type f -delete
find . -name "*.exe" -type f -delete
find . -name "*.x" -type f -delete

# compile the source file and run
g++ main.cpp -o prg.x
./prg.x

# plot the results according to the mode user used
python plot.py