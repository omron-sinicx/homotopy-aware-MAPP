find ./ -regex "./src/.*\.\(c\|h\)pp" | xargs clang-format -style=file -i; find . -name '*~' -exec rm {} \;
