def damerau_levenshtein_dp_optimized(s1: str, s2: str) -> int:
    """
    Optimized dynamic programming implementation of Damerau-Levenshtein distance 
    using a (3 x m+1) matrix.

    Args:
        s1 (str): The first word
        s2 (str): The second word

    Returns:
        int: The Damerau-Levenshtein Distance
    """
    n = len(s1)
    m = len(s2)
    
    # (3 x (m + 1)) matrix to store values
    matrix = [[0] * (m + 1) for _ in range(3)]
    
    # filling the matrix
    for i in range(n + 1):
        for j in range(m + 1):
            # if min(i, j) = 0 return max(i, j)
            if i == 0:
                matrix[i % 3][j] = j
            elif j == 0:
                matrix[i % 3][j] = i

            # otherwise
            else:
                # Extra edit count for substitution depending whether two letters are same
                edit_substitution = 0 if s1[i - 1] == s2[j - 1] else 1
                
                # if it is possible to include transpositon to the possible paths 
                if i > 1 and j > 1 and s1[i - 1] == s2[j - 2] and s1[i - 2] == s2[j - 1]:
                    matrix[i % 3][j] = min(
                        matrix[(i - 1) % 3][j] + 1,                      # Deletion
                        matrix[i % 3][j - 1] + 1,                        # Insertion
                        matrix[(i - 1) % 3][j - 1] + edit_substitution,  # Substitution
                        matrix[(i - 2) % 3][j - 2] + 1                   # Transposition
                    )

                else:
                    matrix[i % 3][j] = min(
                        matrix[(i - 1) % 3][j] + 1,                      # Deletion
                        matrix[i % 3][j - 1] + 1,                        # Insertion
                        matrix[(i - 1) % 3][j - 1] + edit_substitution   # Substitution
                    )
    
    return matrix[-1][-1]


def jaro_similarity(s1: str, s2: str) -> float:
    """
    Computes the Jaro similarity between two strings.

    Args:
        s1 (str): The first string.
        s2 (str): The second string.

    Returns:
        float: Jaro similarity score (0.0 to 1.0)
    """
    n = len(s1)
    m = len(s2)
     
    # complete similar words
    if n == 0 and m == 0:
        return 1.0
    
    match_distance = max(n, m) // 2 - 1

    # matching trackers
    s1_matches = [False] * n
    s2_matches = [False] * m

    # counts
    matches = 0
    transpositions = 0
    
    # Find matches
    for i in range(n):
        # range of positions for matching strings
        start = max(0, i - match_distance)
        end = min(m, i + match_distance + 1)

        # iterate between boundaries of mathcing distance
        for j in range(start, end):
            # if a letter has already been assigned as matched, skip
            if s2_matches[j]:
                continue

            # if letters are differet, skip
            if s1[i] != s2[j]:
                continue
            
            s1_matches[i] = True
            s2_matches[j] = True
            matches += 1

            break
    
    # if there is no matching letters, then no need to check for transpositions.
    if matches == 0:
        return 0.0
    
    # Count transpositions
    k = 0
    # iterate through letters of 1st string
    for i in range(n):
        # if transposition is possible, then it should be in the matching tracker
        if not s1_matches[i]:
            continue
        
        # iterating until matching letter found in the 2nd word
        while not s2_matches[k]:
            k += 1
        
        if s1[i] != s2[k]:
            transpositions += 1
        k += 1

    transpositions //= 2

    jaro_sim = (matches / n + matches / m + (matches - transpositions) / matches) / 3.0

    return jaro_sim


def jaro_winkler_similarity(s1: str, s2: str, scaling_factor: float = 0.1) -> float:
    """
    Computes the Jaro-Winkler similarity between two strings.

    Args:
        s1 (str): The first string.
        s2 (str): The second string.
        scaling_factor (float): The weight for common prefixes (default is 0.1).

    Returns:
        float: Jaro-Winkler similarity score (0.0 to 1.0)
    """
    
    jaro_sim = jaro_similarity(s1, s2)
    prefix_length = 0
    
    # Calculate common prefix length (max 4 characters)
    for i in range(min(len(s1), len(s2), 4)):
        if s1[i] == s2[i]:
            prefix_length += 1
        
        # ending iteration when mismatch of letters has been faced
        else:
            break

    jw = jaro_sim + prefix_length * scaling_factor * (1 - jaro_sim)
    
    return jw





def hamming_distance(s1: str, s2: str) -> int:
    """
    Computes the Hamming distance between two equal-length strings.

    Args:
        s1 (str): The first string.
        s2 (str): The second string.

    Returns:
        int: Hamming distance.
    """
    if len(s1) != len(s2):
        raise ValueError("Strings must be of equal length for Hamming distance.")
    
    distance = 0

    for x, y in zip(s1, s2):
        if x != y:
            distance += 1

    return distance

def hamming_distance_binary(s1: str, s2: str) -> int:
    """
    Calculate the Hamming Distance between two strings.

    Args:
        s1 (str): The first string to compare.
        s2 (str): The second string to compare.

    Returns:
        int: The Hamming Distance between the two strings.
    """
    if len(s1) != len(s2):
        raise ValueError("Strings must be of equal length for Hamming distance.")
    
    distance = 0
    for x, y in zip(s1, s2):
        # Convert each character to its binary form (ASCII/Unicode) and compare
        bin_x = format(ord(x), '08b')
        bin_y = format(ord(y), '08b')
        
        # Count differing bits
        for b1, b2 in zip(bin_x, bin_y):
            if b1 != b2:
                distance += 1
    
    return distance