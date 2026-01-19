wolframscript -f bj2.m 
cat "kxx_shear.txt" | sed 's/List(/kxx_shear = {/' | sed 's/)$/};/'  > kxx_shear.cpp
cat "knx_shear.txt" | sed 's/List(/knx_shear = {/' | sed 's/)$/};/'  > knx_shear.cpp
cat "knn_shear.txt" | sed 's/^/knn_shear = /' | sed 's/$/;/'  > knn_shear.cpp
rm knx_shear.txt kxx_shear.txt knn_shear.txt

cat "kxx_bulk.txt" | sed 's/List(/kxx_bulk = {/' | sed 's/)$/};/'  > kxx_bulk.cpp
cat "knx_bulk.txt" | sed 's/List(/knx_bulk = {/' | sed 's/)$/};/'  > knx_bulk.cpp
cat "knn_bulk.txt" | sed 's/^/knn_bulk = /' | sed 's/$/;/'  > knn_bulk.cpp
rm knx_bulk.txt kxx_bulk.txt knn_bulk.txt

cp kxx_shear.cpp ./VischydroNode_inc.hpp 
echo "\n" >> ./VischydroNode_inc.hpp
cat knx_shear.cpp >> ./VischydroNode_inc.hpp
echo "\n" >> ./VischydroNode_inc.hpp
cat knn_shear.cpp >> ./VischydroNode_inc.hpp
echo "\n" >> ./VischydroNode_inc.hpp
cat kxx_bulk.cpp >> ./VischydroNode_inc.hpp
echo "\n">> ./VischydroNode_inc.hpp
cat knx_bulk.cpp >> ./VischydroNode_inc.hpp
echo "\n">> ./VischydroNode_inc.hpp
cat knn_bulk.cpp >> ./VischydroNode_inc.hpp  
echo "\n">> ./VischydroNode_inc.hpp
rm kxx_shear.cpp knx_shear.cpp knn_shear.cpp kxx_bulk.cpp knx_bulk.cpp knn_bulk.cpp
