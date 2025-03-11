def levenshtein_recursive(s1: str, s2: str) -> int: 
    """
    Recursive implementation of Levenshtein distance.

    Args:
        s1 (str): The first word
        s2 (str): The second word

    Returns:
        int: The Levenshtein Distance
    """

    # Condition return max(i, j) if min(i, j) = 0
    if not s1:
        return len(s2)
    
    if not s2:
        return len(s1)
    
    # otherwise
    # Extra edit count for substitution depending whether two letters are same
    edit_substitution = 0 if s1[0] == s2[0] else 1
    return min(
        levenshtein_recursive(s1[1:], s2) + 1,                          # Deletion
        levenshtein_recursive(s1, s2[1:]) + 1,                          # Insertion
        levenshtein_recursive(s1[1:], s2[1:]) + edit_substitution       # Substitution
    )


def levenshtein_dp(s1: str, s2: str) -> tuple[int, list[list[int]]]:
    """
    Dynamic programming implementation of Levenshtein distance.

    Args:
        s1 (str): The first word
        s2 (str): The second word

    Returns:
        tuple[int, list[list[int]]]: 
            - The Levenshtein Distance (int).
            - The Matrix of the Levenshtein Distance (list of lists of integers).
    """

    n = len(s1)
    m = len(s2)

    # (n + 1, m + 1) matrix to store values
    matrix = [[0] * (m + 1) for _ in range(n + 1)]

    # filling matrix
    for i in range(n + 1):
        for j in range(m + 1):
            # if min(i, j) = 0 return max(i, j)
            if i == 0:
                matrix[i][j] = j

            elif j == 0:
                matrix[i][j] = i

            # otherwise
            else:
                # Extra edit count for substitution depending whether two letters are same
                edit_substitution = 0 if s1[i - 1] == s2[j - 1] else 1

                matrix[i][j] = min(
                    matrix[i - 1][j] + 1,                        # Deletion
                    matrix[i][j - 1] + 1,                        # Insertion
                    matrix[i - 1][j - 1] + edit_substitution     # Substitution
                )
    
    return matrix[-1][-1], matrix

def levenshtein_dp_optimized(s1: str, s2: str) -> int:
    """
    Optimized dynamic programming implementation of Levenshtein distance 
    by using a (2 x m+1) matrix instead of (n+1 x m+1) to store values.

    Args:
        s1 (str): The first word
        s2 (str): The second word

    Returns:
        int: The Levenshtein Distance
    """

    n = len(s1)
    m = len(s2)

    # (2 x (m + 1)) matrix to store values
    matrix = [[0] * (m + 1) for _ in range(2)]

    # filling the matrix
    for i in range(n + 1):
        for j in range(m + 1):
            # if min(i, j) = 0 return max(i, j)
            if i == 0:
                matrix[i % 2][j] = j
            elif j == 0:
                matrix[i % 2][j] = i

            # otherwise
            else:
                # Extra edit count for substitution depending whether two letters are same
                edit_substitution = 0 if s1[i - 1] == s2[j - 1] else 1

                matrix[i % 2][j] = min(
                    matrix[(i - 1) % 2][j] + 1,                      # Deletion
                    matrix[i % 2][j - 1] + 1,                        # Insertion
                    matrix[(i - 1) % 2][j - 1] + edit_substitution   # Substitution
                )

    return matrix[-1][-1]