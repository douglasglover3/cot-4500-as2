import numpy as np
np.set_printoptions(precision=7, suppress=True, linewidth=100)

# Question 1

x = [3.6, 3.8, 3.9]
fx = [1.675, 1.436, 1.318]
size = 3

desiredX = 3.7
sum = 0
for k in range(size):
    product = 1
    for i in range(size):
        if(i != k):
            product *= (desiredX - x[i]) / (x[k] - x[i])
    sum += product * fx[k]

print(sum, end="\n\n")

# Question 2

# point setup
x_points = [7.2, 7.4, 7.5, 7.6]
y_points = [23.5492, 25.3913, 26.8224, 27.4589]

def divided_difference_table(x_points, y_points):
    # set up the matrix
    size: int = np.size(x_points)
    matrix: np.array = np.zeros((size, size))
    # fill the matrix
    for index, row in enumerate(matrix):
        row[0] =  y_points[index]
    # populate the matrix (end points are based on matrix size and max operations we're using)
    for i in range(1, size):
        for j in range(1, size):
            # the numerator are the immediate left and diagonal left indices...
            numerator = matrix[i][j - 1] - matrix[i - 1][j - 1]
            # the denominator is the X-SPAN...
            denominator = x_points[i] - x_points[i - j]
            operation = numerator / denominator
            # cut it off to view it more simpler
            matrix[i][j] = operation
    return matrix

divided_table = divided_difference_table(x_points, y_points)
print([divided_table[1][1], divided_table[2][2], divided_table[3][3]], end="\n\n")


# Question 3

def get_approximate_result(matrix, x_points, value):
    # p0 is always y0 and we use a reoccuring x to avoid having to recalculate x 
    reoccuring_x_span = 1
    reoccuring_px_result = matrix[0][0]
    
    # we only need the diagonals...and that starts at the first row...
    for index in range(1, np.size(x_points)):
        polynomial_coefficient = matrix[index][index]
        # we use the previous index for x_points....
        reoccuring_x_span *= (value - x_points[index - 1])
        
        # get a_of_x * the x_span
        mult_operation = polynomial_coefficient * reoccuring_x_span
        # add the reoccuring px result
        reoccuring_px_result += mult_operation
    
    # final result
    return reoccuring_px_result

# find approximation
approximating_x = 7.3
final_approximation = get_approximate_result(divided_table, x_points, approximating_x)

print(final_approximation, end="\n\n")


# Question 4

def divided_difference_hermite_table(x_points, y_points, prime_points):
    # set up the matrix
    size: int = np.size(x_points * 2)
    matrix: np.array = np.zeros((size, size))
    # fill the matrix
    for index, row in enumerate(matrix):
        row[0] =  x_points[int(index / 2)]
        row[1] =  y_points[int(index / 2)]
        if index % 2 != 0:
            row[2] = prime_points[int(index / 2)]
    # populate the matrix (end points are based on matrix size and max operations we're using)
    for i in range(1, size):
        for j in range(1, i + 1):
            if(j >= size - 1):
                continue
            # the numerator are the immediate left and diagonal left indices...
            numerator = matrix[i][j] - matrix[i - 1][j]
            # the denominator is the X-SPAN...
            denominator = x_points[int(i / 2)] - x_points[int((i - j) / 2)]
            if denominator > 0:
                operation = numerator / denominator
                # cut it off to view it more simpler
                matrix[i][j + 1] = operation
            
    return matrix

# point setup
x = [3.6, 3.8, 3.9]
y = [1.675, 1.436, 1.318]
yprime = [-1.195, -1.188, -1.182]
divided_table = divided_difference_hermite_table(x, y, yprime)
print(divided_table, end="\n\n")


# Question 5
x = [2, 5, 8, 10]
y = [3, 5, 7, 9]

def cubicSplineInterpolation(x, a):
    n = np.size(x) - 1

    h = np.zeros(n)
    for i in range(n):
        h[i] = x[i + 1] - x[i]

    matrixA = np.zeros((n + 1, n + 1))
    matrixA[0][0] = 1
    matrixA[n][n] = 1
    for i in range(1, n):
        matrixA[i][i - 1] = h[i - 1]
        matrixA[i][i] = 2 * (h[i - 1] + h[i])
        matrixA[i][i + 1] = h[i]
    print(matrixA, end="\n\n")

    n = np.size(a)
    b = np.zeros(n)
    for i in range(2, n):
        b[i - 1] = ((3 / h[i - 1]) * (a[i] - a[i - 1])) - ((3 / h[i - 2]) * (a[i - 1] - a[i - 2]))
    print(b, end="\n\n")

    l = np.ones(n)
    u = np.zeros(n)
    z = np.zeros(n)
    c = np.zeros(n)
    for i in range(1, n - 1):
        l[i] = 2 * (x[i + 1] - x[i - 1]) - h[i - 1] * u[i - 1]
        u[i] = h[i] / l[i]
        z[i] = (b[i] - h[i - 1] * z[i - 1]) / l[i]
    
    for j in np.flip(range(1, n - 1)):
        c[j] = z[j] - u[j] * c[j + 1]
    print(c, end="\n\n")

cubicSplineInterpolation(x, y)