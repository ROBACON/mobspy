import re


# Compare results with expected file
def compare_model(comp_results, file_name):
    def normalize(text):
        """Normalize text by removing extra spaces, tabs, and Unicode variations."""
        text = text.strip()  # Remove leading/trailing spaces
        text = re.sub(
            r"\s+", " ", text
        )  # Replace multiple spaces/tabs/newlines with a single space
        return text

    with open(file_name, "r", encoding="utf-8") as file:
        for r, line in zip(comp_results.split("\n"), file.readlines()):
            line = normalize(line)
            r = normalize(r)

            if r != line:
                print("file: " + line)
                print("test: " + r)
                return False
    return True


def compare_model_ignore_order(comp_results, file_name):
    with open(file_name, "r") as file:
        file_lines = {line.strip() for line in file.readlines() if line.strip()}
        result_lines = {
            line.strip() for line in comp_results.split("\n") if line.strip()
        }

        # Check if both sets of lines are equal
        if file_lines != result_lines:
            missing_in_file = result_lines - file_lines
            missing_in_results = file_lines - result_lines

            if missing_in_file:
                print("Lines in results but not in file:")
                for line in missing_in_file:
                    print(f"test: {line}")

            if missing_in_results:
                print("Lines in file but not in results:")
                for line in missing_in_results:
                    print(f"file: {line}")

            return False
    return True
