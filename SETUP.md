
# Set up MacOS example

In case you need to get pyMultiNest to work on a Mac computer.

## get Multinest

Pick a `/path/to/MultiNest`, then ...
```
git clone https://github.com/JohannesBuchner/MultiNest.git
cd MultiNest/build/
cmake ..
make
sudo make install
```

## install pyMultinest

```
python -m pip install --force-reinstall pymultinest
MN_DIR=/path/to/MultiNest/ python -m pip install --force-reinstall pymultinest
```

## run pyMultinest

Before starting Python, do
```
export DYLD_LIBRARY_PATH=/path/to/MultiNest/lib:$DYLD_LIBRARY_PATH
```
This command needs to be run every time a new terminal is opened.

