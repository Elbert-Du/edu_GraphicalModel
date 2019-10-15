rm *.so
rm privBayes.cpp
rm -r build

python3 setup.py build_ext --inplace 
