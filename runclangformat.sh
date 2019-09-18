find ./test -name *.h -type f -exec sed -i 's/\t/    /;s/[[:space:]]*$//' {} \;
find ./test -name *.cc -type f -exec sed -i 's/\t/    /;s/[[:space:]]*$//' {} \;
find ./test -name *.c -type f -exec sed -i 's/\t/    /;s/[[:space:]]*$//' {} \;
find ./Source -name *.c -type f -exec sed -i 's/\t/    /;s/[[:space:]]*$//' {} \;
clang-format-5.0 -i -style=file ./test/*.h ./test/*.cc ./test/e2e_test/*.h ./test/e2e_test/*.cc ./test/api_test/*.h ./test/api_test/*.cc