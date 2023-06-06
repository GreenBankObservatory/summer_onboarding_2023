import argparse

parser = argparse.ArgumentParser()
parser.add_argument("string")
args = parser.parse_args()
input = str(args.string)


def is_palindrome(string):
    """Determines whether input `string` is a palindrome"""

    front = 0
    end = len(string) - 1
    while front < end:
        if string[front].upper() != string[end].upper():
            return False
        front += 1
        end -= 1
    return True

print(is_palindrome(input))
if is_palindrome(input):
    raise SystemExit(0)
else:
    raise SystemExit(1)

tests = [
    ("askdfjs", False),
    ("racecar", True),
    ("RaceCar", True),
    ("palindrome", False),
    ("1234321", True),
    ("⭐📡👀", False),
    ("🌿🌼🌿", True),
]


def test_is_palindrome():
    failures = []
    for string, expected_result in tests:
        actual_result = is_palindrome(string)
        if actual_result != expected_result:
            failures.append(
                f"Expected is_palindrome({string!r}) to return "
                f"{expected_result}; got {actual_result}"
            )

    if failures:
        fstr = "\n".join(failures)
        raise AssertionError(f"Your function failed with the following inputs:\n{fstr}")
    else:
        print("Success!")


if __name__ == "__main__":
    test_is_palindrome()