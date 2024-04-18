test_set = {1, 2, 3, 4}
print(test_set)

for number in list(test_set):
    print(number)
    test_set.discard(number + 1)

print(test_set)


test_set = {1, 2, 3, 4}
print(test_set)

while test_set:
    number = test_set.pop()
    print(number)
    test_set.discard(number + 1)

print(test_set)

test_set = {1, 2, 3, 4}
print(test_set)
test_set_copy = set(test_set)
test_set_copy.discard(2)
print(test_set)
