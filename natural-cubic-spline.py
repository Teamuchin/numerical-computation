import numpy as np
import matplotlib.pyplot as plt


def natural_cubic_spline(x_points, y_points):

    #Finds the abcd values in an array for each segment of the spline.

    # I'm using matrix method here to derive sufficient functions
    n = len(x_points)

    # Step 1: Calculate delta_x
    # delta_x = x_{i+1} - x_i
    delta_x = [x_points[i + 1] - x_points[i] for i in range(n - 1)]

    # Step 2: Calculate Delta_y
    # Delta_y = (y_{i+1} - y_i) / delta_i
    delta_y = [(y_points[i + 1] - y_points[i]) / delta_x[i] for i in range(n - 1)]

    # Step 3: Set up the system of equations for c_i coefficients
    # There will be n c_i values (c_0 to c_{n-1})
    # The matrix A will be (n x n) and vector B will be (n x 1)

    A = np.zeros((n, n))
    B = np.zeros(n)

    # First row: c_0 = 0 for natural spline
    A[0, 0] = 1
    B[0] = 0

    # Last row: c_{n-1} = 0 for natural spline
    A[n - 1, n - 1] = 1
    B[n - 1] = 0

    # Fill in the intermediate rows (1 to n-2)
    for i in range(1, n - 1):
        A[i, i - 1] = delta_x[i - 1]  # delta_n at nth column
        A[i, i] = 2 * (delta_x[i - 1] + delta_x[i])  # 2 * (delta_n + delta_{n+1}) on n+1th column
        A[i, i + 1] = delta_x[i]  # delta_{n+1} at n+2th column
        B[i] = 3 * (delta_y[i] - delta_y[i - 1])  # 3 * ( (Delta_{n+1}/delta_{n+1}) - (Delta_n/delta_n) )

    c_coeffs = np.linalg.solve(A, B)


    # Calculate a, b, d coefficients for each segment
    a_coeffs = [y_points[i] for i in range(n - 1)]  # S_i(x_i) = y_i, so a_i = y_i
    b_coeffs = []
    d_coeffs = []

    for i in range(n - 1):
        # b_i = (y_{i+1} - y_i) / delta_i - delta_i * (2*c_i + c_{i+1}) / 3
        b_val = delta_y[i] - (delta_x[i] * (2 * c_coeffs[i] + c_coeffs[i + 1])) / 3
        b_coeffs.append(b_val)

        # d_i = (c_{i+1} - c_i) / (3 * delta_i)
        d_val = (c_coeffs[i + 1] - c_coeffs[i]) / (3 * delta_x[i])
        d_coeffs.append(d_val)

    return a_coeffs, b_coeffs, c_coeffs.tolist()[
                               :-1], d_coeffs  # c_coeffs has n elements, we need n-1 for S_i(x-x_i)^2 term


def evaluate_natural_cubic_spline(x, x_points, a_coeffs, b_coeffs, c_coeffs, d_coeffs):

    #Evaluates the natural cubic spline at a given x-value.


    # Find the correct segment
    # The segments are [x_i, x_{i+1}]
    segment_index = -1
    for i in range(len(x_points) - 1):
        if x_points[i] <= x <= x_points[i + 1]:
            segment_index = i
            break

    # Handle the case where x is exactly the last point
    if x == x_points[-1] and segment_index == -1 and len(x_points) > 1:
        segment_index = len(x_points) - 2  # Use the last segment

    if segment_index == -1:
        # x is outside the range of x_points or no valid segment found
        # For interpolation, we generally don't extrapolate.
        return None

    # Get coefficients for the chosen segment
    a = a_coeffs[segment_index]
    b = b_coeffs[segment_index]
    c = c_coeffs[segment_index]
    d = d_coeffs[segment_index]

    # Calculate (x - x_i) for the segment
    h = x - x_points[segment_index]

    # Evaluate the polynomial S_i(x) = a_i + b_i(x-x_i) + c_i(x-x_i)^2 + d_i(x-x_i)^3
    y = a + b * h + c * (h ** 2) + d * (h ** 3)
    return y

#face circle cords which i have get help from desmos to obtain coords of the points i used desmos on finding cords for other circular shapes too
x_data = [-2,-1.975,-1.875,-1.75,-1.5,-1.25,-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1,1.25,1.5,1.75,1.875,1.975,2]
y_data = [0,0.31524,0.69597, 0.96825, 1.32288, 1.56125, 1.73205, 1.85405, 1.93649, 1.98431, 2, 1.98431, 1.93649, 1.85405, 1.73205, 1.56125, 1.32288, 0.96825,0.69597,0.31524, 0]
x_data_1 = [-2,-1.975,-1.875,-1.75,-1.5,-1.25,-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1,1.25,1.5,1.75,1.875,1.975,2]
y_data_1 = [0,-0.31524,-0.69597, -0.96825, -1.32288, -1.56125, -1.73205, -1.85405, -1.93649, -1.98431, -2, -1.98431, -1.93649, -1.85405, -1.73205, -1.56125, -1.32288, -0.96825,-0.69597,-0.31524, 0]
#other face cords
x_data_2 = [-1,-0.75,-0.5]
y_data_2 = [1,1.2,1]
x_data_3 = [0.5,0.75,1]
y_data_3 = [1,1.2,1]
x_data_4 = [-1,0,1]
y_data_4 = [-0.5,-1,-0.5]
x_data_5 = [-1.25,-1.125,-1]
y_data_5 = [0.25,0.375,0.25]
x_data_6 = [1,1.125,1.25]
y_data_6 = [0.25,0.375,0.25]

a, b, c, d = natural_cubic_spline(x_data, y_data)
a1,b1,c1,d1 = natural_cubic_spline(x_data_1, y_data_1)
a2,b2,c2,d2 = natural_cubic_spline(x_data_2, y_data_2)
a3,b3,c3,d3 = natural_cubic_spline(x_data_3, y_data_3)
a4,b4,c4,d4 = natural_cubic_spline(x_data_4, y_data_4)
a5,b5,c5,d5 = natural_cubic_spline(x_data_5, y_data_5)
a6,b6,c6,d6 = natural_cubic_spline(x_data_6, y_data_6)

x_values = np.linspace(min(x_data), max(x_data), 50)
x_values_1 = np.linspace(min(x_data_1), max(x_data_1), 50)
x_values_2 = np.linspace(min(x_data_2), max(x_data_2), 50)
x_values_3 = np.linspace(min(x_data_3), max(x_data_3), 50)
x_values_4 = np.linspace(min(x_data_4), max(x_data_4), 50)
x_values_5 = np.linspace(min(x_data_5), max(x_data_5), 50)
x_values_6 = np.linspace(min(x_data_6), max(x_data_6), 50)

interpolated_y_values = []
interpolated_y_values_1 = []
interpolated_y_values_2 = []
interpolated_y_values_3 = []
interpolated_y_values_4 = []
interpolated_y_values_5 = []
interpolated_y_values_6 = []

for x_val in x_values:
    y_val = evaluate_natural_cubic_spline(x_val, x_data, a, b, c, d)
    interpolated_y_values.append(y_val)

for x_val_1 in x_values_1:
    y_val_1 = evaluate_natural_cubic_spline(x_val_1, x_data_1, a1, b1, c1, d1)
    interpolated_y_values_1.append(y_val_1)

for x_val_2 in x_values_2:
    y_val_2 = evaluate_natural_cubic_spline(x_val_2, x_data_2, a2, b2, c2, d2)
    interpolated_y_values_2.append(y_val_2)

for x_val_3 in x_values_3:
    y_val_3 = evaluate_natural_cubic_spline(x_val_3, x_data_3, a3, b3, c3, d3)
    interpolated_y_values_3.append(y_val_3)

for x_val_4 in x_values_4:
    y_val_4 = evaluate_natural_cubic_spline(x_val_4, x_data_4, a4, b4, c4, d4)
    interpolated_y_values_4.append(y_val_4)

for x_val_5 in x_values_5:
    y_val_5 = evaluate_natural_cubic_spline(x_val_5, x_data_5, a5, b5, c5, d5)
    interpolated_y_values_5.append(y_val_5)

for x_val_6 in x_values_6:
    y_val_6 = evaluate_natural_cubic_spline(x_val_6, x_data_6, a6, b6, c6, d6)
    interpolated_y_values_6.append(y_val_6)


plt.figure(figsize=(10,10))
plt.plot(x_data, y_data, 'o', color='blue', markersize=8)
plt.plot(x_data_1, y_data_1, 'o', color='red', markersize=8)
plt.plot(x_data_2, y_data_2, 'o', color='green', markersize=8)
plt.plot(x_data_3, y_data_3, 'o', color='purple', markersize=8)
plt.plot(x_data_4, y_data_4, 'o', color='orange', markersize=8)
plt.plot(x_data_5, y_data_5, 'o', color='black', markersize=8)
plt.plot(x_data_6, y_data_6, 'o', color='yellow', markersize=8)


plt.plot(x_values, interpolated_y_values, '-', color='blue', linewidth=2)
plt.plot(x_values_1, interpolated_y_values_1, '-', color='red', linewidth=2)
plt.plot(x_values_2, interpolated_y_values_2, '-', color='green', linewidth=2)
plt.plot(x_values_3, interpolated_y_values_3, '-', color='purple', linewidth=2)
plt.plot(x_values_4, interpolated_y_values_4, '-', color='orange', linewidth=2)
plt.plot(x_values_5, interpolated_y_values_5, '-', color='black', linewidth=2)
plt.plot(x_values_6, interpolated_y_values_6, '-', color='yellow', linewidth=2)

plt.title('Smile With Blush')
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.grid(True)
plt.axis('equal')
plt.show()




"""
Next face
"""


#other face cords
x_data_2 = [-1,-0.5]
y_data_2 = [1.25,1]
x_data_3 = [0.5,1]
y_data_3 = [1,1.25]
x_data_4 = [-1,0,1]
y_data_4 = [-1,-0.5,-1]

a2,b2,c2,d2 = natural_cubic_spline(x_data_2, y_data_2)
a3,b3,c3,d3 = natural_cubic_spline(x_data_3, y_data_3)
a4,b4,c4,d4 = natural_cubic_spline(x_data_4, y_data_4)

x_values_2 = np.linspace(min(x_data_2), max(x_data_2), 50)
x_values_3 = np.linspace(min(x_data_3), max(x_data_3), 50)
x_values_4 = np.linspace(min(x_data_4), max(x_data_4), 50)

interpolated_y_values_2 = []
interpolated_y_values_3 = []
interpolated_y_values_4 = []


for x_val_2 in x_values_2:
    y_val_2 = evaluate_natural_cubic_spline(x_val_2, x_data_2, a2, b2, c2, d2)
    interpolated_y_values_2.append(y_val_2)

for x_val_3 in x_values_3:
    y_val_3 = evaluate_natural_cubic_spline(x_val_3, x_data_3, a3, b3, c3, d3)
    interpolated_y_values_3.append(y_val_3)

for x_val_4 in x_values_4:
    y_val_4 = evaluate_natural_cubic_spline(x_val_4, x_data_4, a4, b4, c4, d4)
    interpolated_y_values_4.append(y_val_4)



plt.figure(figsize=(10,10))

plt.plot(x_data, y_data, 'o', color='blue', markersize=8)
plt.plot(x_data_1, y_data_1, 'o', color='red', markersize=8)
plt.plot(x_data_2, y_data_2, 'o', color='green', markersize=8)
plt.plot(x_data_3, y_data_3, 'o', color='purple', markersize=8)
plt.plot(x_data_4, y_data_4, 'o', color='orange', markersize=8)

plt.plot(x_values, interpolated_y_values, '-', color='blue', linewidth=2)
plt.plot(x_values_1, interpolated_y_values_1, '-', color='red', linewidth=2)
plt.plot(x_values_2, interpolated_y_values_2, '-', color='green', linewidth=2)
plt.plot(x_values_3, interpolated_y_values_3, '-', color='purple', linewidth=2)
plt.plot(x_values_4, interpolated_y_values_4, '-', color='orange', linewidth=2)

plt.title('Angry')
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.grid(True)
plt.axis('equal')
plt.show()


"""
Next face
"""


x_data_2 = [-1,-0.75,-0.5]
y_data_2 = [1,1.2,1]
x_data_3 = [0.5,0.75,1]
y_data_3 = [1,1.2,1]
x_data_4 = [-1,0,1]
y_data_4 = [-0.5,-1,-0.5]
x_data_5 = [-1.5,-1.25,-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1,1.25,1.5]
y_data_5 = [2.01118, 2.276542, 2.372753, 2.433049, 2.471419, 2.49301, 2.5,2.49301, 2.471419, 2.433049, 2.372753, 2.276542, 2.01118]
x_data_6 = [-1.5,-1.25,-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1,1.25,1.5]
y_data_6 = [1.98882, 1.723458, 1.627247, 1.566951, 1.528581, 1.50699, 1.5,1.50699, 1.528581, 1.566951, 1.627247, 1.723458, 1.98882]


a2,b2,c2,d2 = natural_cubic_spline(x_data_2, y_data_2)
a3,b3,c3,d3 = natural_cubic_spline(x_data_3, y_data_3)
a4,b4,c4,d4 = natural_cubic_spline(x_data_4, y_data_4)
a5,b5,c5,d5 = natural_cubic_spline(x_data_5, y_data_5)
a6,b6,c6,d6 = natural_cubic_spline(x_data_6, y_data_6)

x_values_2 = np.linspace(min(x_data_2), max(x_data_2), 50)
x_values_3 = np.linspace(min(x_data_3), max(x_data_3), 50)
x_values_4 = np.linspace(min(x_data_4), max(x_data_4), 50)
x_values_5 = np.linspace(min(x_data_5), max(x_data_5), 50)
x_values_6 = np.linspace(min(x_data_6), max(x_data_6), 50)

interpolated_y_values_2 = []
interpolated_y_values_3 = []
interpolated_y_values_4 = []
interpolated_y_values_5 = []
interpolated_y_values_6 = []


for x_val_2 in x_values_2:
    y_val_2 = evaluate_natural_cubic_spline(x_val_2, x_data_2, a2, b2, c2, d2)
    interpolated_y_values_2.append(y_val_2)

for x_val_3 in x_values_3:
    y_val_3 = evaluate_natural_cubic_spline(x_val_3, x_data_3, a3, b3, c3, d3)
    interpolated_y_values_3.append(y_val_3)

for x_val_4 in x_values_4:
    y_val_4 = evaluate_natural_cubic_spline(x_val_4, x_data_4, a4, b4, c4, d4)
    interpolated_y_values_4.append(y_val_4)

for x_val_5 in x_values_5:
    y_val_5 = evaluate_natural_cubic_spline(x_val_5, x_data_5, a5, b5, c5, d5)
    interpolated_y_values_5.append(y_val_5)

for x_val_6 in x_values_6:
    y_val_6 = evaluate_natural_cubic_spline(x_val_6, x_data_6, a6, b6, c6, d6)
    interpolated_y_values_6.append(y_val_6)


plt.figure(figsize=(10,10))
plt.plot(x_data, y_data, 'o', color='blue',  markersize=8)
plt.plot(x_data_1, y_data_1, 'o', color='red',  markersize=8)
plt.plot(x_data_2, y_data_2, 'o', color='green',  markersize=8)
plt.plot(x_data_3, y_data_3, 'o', color='purple',  markersize=8)
plt.plot(x_data_4, y_data_4, 'o', color='orange',  markersize=8)
plt.plot(x_data_5, y_data_5, 'o', color='black', markersize=8)
plt.plot(x_data_6, y_data_6, 'o', color='yellow',  markersize=8)



plt.plot(x_values, interpolated_y_values, '-', color='blue', linewidth=2)
plt.plot(x_values_1, interpolated_y_values_1, '-', color='red', linewidth=2)
plt.plot(x_values_2, interpolated_y_values_2, '-', color='green', linewidth=2)
plt.plot(x_values_3, interpolated_y_values_3, '-', color='purple', linewidth=2)
plt.plot(x_values_4, interpolated_y_values_4, '-', color='orange', linewidth=2)
plt.plot(x_values_5, interpolated_y_values_5, '-', color='black', linewidth=2)
plt.plot(x_values_6, interpolated_y_values_6, '-', color='yellow', linewidth=2)

plt.title('Smile With Halo')
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.grid(True)
plt.axis('equal')
plt.show()


"""
Letters
"""


x_data = [-2,2]
y_data = [2,2]
x_data_1 = [0,0.0001]
y_data_1 = [2,-2]


a,b,c,d = natural_cubic_spline(x_data, y_data)
a1,b1,c1,d1 = natural_cubic_spline(x_data_1, y_data_1)

x_values = np.linspace(min(x_data), max(x_data), 50)
x_values_1 = np.linspace(min(x_data_1), max(x_data_1), 50)


interpolated_y_values = []
interpolated_y_values_1 = []


for x_val in x_values:
    y_val = evaluate_natural_cubic_spline(x_val, x_data, a, b, c, d)
    interpolated_y_values.append(y_val)

for x_val_1 in x_values_1:
    y_val_1 = evaluate_natural_cubic_spline(x_val_1, x_data_1, a1, b1, c1, d1)
    interpolated_y_values_1.append(y_val_1)


plt.figure(figsize=(10,10))
plt.plot(x_data, y_data, 'o', color='blue',  markersize=8)
plt.plot(x_data_1, y_data_1, 'o', color='red',  markersize=8)


plt.plot(x_values, interpolated_y_values, '-', color='blue', linewidth=2)
plt.plot(x_values_1, interpolated_y_values_1, '-', color='red', linewidth=2)

plt.title('Letter T')
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.grid(True)
plt.axis('equal')
plt.show()


"""
Letters
"""

x_data = [-2,-1.975,-1.875,-1.75,-1.5,-1.25,-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1,1.25,1.5,1.75]
y_data = [0,0.31524,0.69597, 0.96825, 1.32288, 1.56125, 1.73205, 1.85405, 1.93649, 1.98431, 2, 1.98431, 1.93649, 1.85405, 1.73205, 1.56125, 1.32288, 0.96825]
x_data_1 = [-2,-1.975,-1.875,-1.75,-1.5,-1.25,-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1,1.25,1.5,1.75]
y_data_1 = [0,-0.31524,-0.69597, -0.96825, -1.32288, -1.56125, -1.73205, -1.85405, -1.93649, -1.98431, -2, -1.98431, -1.93649, -1.85405, -1.73205, -1.56125, -1.32288, -0.96825]
x_data_2 = [0,0.0001]
y_data_2 = [-2,-2.2]

a,b,c,d = natural_cubic_spline(x_data, y_data)
a1,b1,c1,d1 = natural_cubic_spline(x_data_1, y_data_1)
a2,b2,c2,d2 = natural_cubic_spline(x_data_2, y_data_2)

x_values = np.linspace(min(x_data), max(x_data), 50)
x_values_1 = np.linspace(min(x_data_1), max(x_data_1), 50)
x_values_2 = np.linspace(min(x_data_2), max(x_data_2), 50)


interpolated_y_values = []
interpolated_y_values_1 = []
interpolated_y_values_2 = []


for x_val in x_values:
    y_val = evaluate_natural_cubic_spline(x_val, x_data, a, b, c, d)
    interpolated_y_values.append(y_val)

for x_val_1 in x_values_1:
    y_val_1 = evaluate_natural_cubic_spline(x_val_1, x_data_1, a1, b1, c1, d1)
    interpolated_y_values_1.append(y_val_1)

for x_val_2 in x_values_2:
    y_val_2 = evaluate_natural_cubic_spline(x_val_2, x_data_2, a2, b2, c2, d2)
    interpolated_y_values_2.append(y_val_2)


plt.figure(figsize=(10,10))
plt.plot(x_data, y_data, 'o', color='blue',  markersize=8)
plt.plot(x_data_1, y_data_1, 'o', color='red',  markersize=8)
plt.plot(x_data_2, y_data_2, 'o', color='green',  markersize=8)


plt.plot(x_values, interpolated_y_values, '-', color='blue', linewidth=2)
plt.plot(x_values_1, interpolated_y_values_1, '-', color='red', linewidth=2)
plt.plot(x_values_2, interpolated_y_values_2, '-', color='green', linewidth=2)

plt.title('Letter Ç')
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.grid(True)
plt.axis('equal')
plt.show()


"""
Letters
"""




x_data = [-0.5,0.5]
y_data = [2,2]
x_data_1 = [-0.5,0.5]
y_data_1 = [-2,-2]
x_data_2 = [0,0.0001]
y_data_2 = [2,-2]
x_data_3 = [-0.05,0.05]
y_data_3 = [2.2,2.1]

a,b,c,d = natural_cubic_spline(x_data, y_data)
a1,b1,c1,d1 = natural_cubic_spline(x_data_1, y_data_1)
a2,b2,c2,d2 = natural_cubic_spline(x_data_2, y_data_2)
a3,b3,c3,d3 = natural_cubic_spline(x_data_3, y_data_3)

x_values = np.linspace(min(x_data), max(x_data), 50)
x_values_1 = np.linspace(min(x_data_1), max(x_data_1), 50)
x_values_2 = np.linspace(min(x_data_2), max(x_data_2), 50)
x_values_3 = np.linspace(min(x_data_3), max(x_data_3), 50)


interpolated_y_values = []
interpolated_y_values_1 = []
interpolated_y_values_2 = []
interpolated_y_values_3 = []


for x_val in x_values:
    y_val = evaluate_natural_cubic_spline(x_val, x_data, a, b, c, d)
    interpolated_y_values.append(y_val)

for x_val_1 in x_values_1:
    y_val_1 = evaluate_natural_cubic_spline(x_val_1, x_data_1, a1, b1, c1, d1)
    interpolated_y_values_1.append(y_val_1)

for x_val_2 in x_values_2:
    y_val_2 = evaluate_natural_cubic_spline(x_val_2, x_data_2, a2, b2, c2, d2)
    interpolated_y_values_2.append(y_val_2)

for x_val_3 in x_values_3:
    y_val_3 = evaluate_natural_cubic_spline(x_val_3, x_data_3, a3, b3, c3, d3)
    interpolated_y_values_3.append(y_val_3)


plt.figure(figsize=(10,10))
plt.plot(x_data, y_data, 'o', color='blue',  markersize=8)
plt.plot(x_data_1, y_data_1, 'o', color='red',  markersize=8)
plt.plot(x_data_2, y_data_2, 'o', color='green',  markersize=8)
plt.plot(x_data_3, y_data_3, 'o', color='purple',  markersize=8)


plt.plot(x_values, interpolated_y_values, '-', color='blue', linewidth=2)
plt.plot(x_values_1, interpolated_y_values_1, '-', color='red', linewidth=2)
plt.plot(x_values_2, interpolated_y_values_2, '-', color='green', linewidth=2)
plt.plot(x_values_3, interpolated_y_values_3, '-', color='purple', linewidth=2)

plt.title('Letter İ')
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.grid(True)
plt.axis('equal')
plt.show()


"""
Letters
"""



x_data = [-2,-1.999]
y_data = [-2,2]
x_data_1 = [-1.999,0]
y_data_1 = [2,-2]
x_data_2 = [0,1.999]
y_data_2 = [-2,2]
x_data_3 = [1.999,2]
y_data_3 = [2,-2]

a,b,c,d = natural_cubic_spline(x_data, y_data)
a1,b1,c1,d1 = natural_cubic_spline(x_data_1, y_data_1)
a2,b2,c2,d2 = natural_cubic_spline(x_data_2, y_data_2)
a3,b3,c3,d3 = natural_cubic_spline(x_data_3, y_data_3)

x_values = np.linspace(min(x_data), max(x_data), 50)
x_values_1 = np.linspace(min(x_data_1), max(x_data_1), 50)
x_values_2 = np.linspace(min(x_data_2), max(x_data_2), 50)
x_values_3 = np.linspace(min(x_data_3), max(x_data_3), 50)


interpolated_y_values = []
interpolated_y_values_1 = []
interpolated_y_values_2 = []
interpolated_y_values_3 = []


for x_val in x_values:
    y_val = evaluate_natural_cubic_spline(x_val, x_data, a, b, c, d)
    interpolated_y_values.append(y_val)

for x_val_1 in x_values_1:
    y_val_1 = evaluate_natural_cubic_spline(x_val_1, x_data_1, a1, b1, c1, d1)
    interpolated_y_values_1.append(y_val_1)

for x_val_2 in x_values_2:
    y_val_2 = evaluate_natural_cubic_spline(x_val_2, x_data_2, a2, b2, c2, d2)
    interpolated_y_values_2.append(y_val_2)

for x_val_3 in x_values_3:
    y_val_3 = evaluate_natural_cubic_spline(x_val_3, x_data_3, a3, b3, c3, d3)
    interpolated_y_values_3.append(y_val_3)


plt.figure(figsize=(10,10))
plt.plot(x_data, y_data, 'o', color='blue',  markersize=8)
plt.plot(x_data_1, y_data_1, 'o', color='red',  markersize=8)
plt.plot(x_data_2, y_data_2, 'o', color='green',  markersize=8)
plt.plot(x_data_3, y_data_3, 'o', color='purple',  markersize=8)


plt.plot(x_values, interpolated_y_values, '-', color='blue', linewidth=2)
plt.plot(x_values_1, interpolated_y_values_1, '-', color='red', linewidth=2)
plt.plot(x_values_2, interpolated_y_values_2, '-', color='green', linewidth=2)
plt.plot(x_values_3, interpolated_y_values_3, '-', color='purple', linewidth=2)

plt.title('Letter M')
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.grid(True)
plt.axis('equal')
plt.show()


"""
Letters
"""




x_data = [-1.5,-1.499]
y_data = [-2,2]
x_data_1 = [-1.499,1.499]
y_data_1 = [2,-2]
x_data_2 = [1.499,1.5]
y_data_2 = [-2,2]

a,b,c,d = natural_cubic_spline(x_data, y_data)
a1,b1,c1,d1 = natural_cubic_spline(x_data_1, y_data_1)
a2,b2,c2,d2 = natural_cubic_spline(x_data_2, y_data_2)

x_values = np.linspace(min(x_data), max(x_data), 50)
x_values_1 = np.linspace(min(x_data_1), max(x_data_1), 50)
x_values_2 = np.linspace(min(x_data_2), max(x_data_2), 50)

interpolated_y_values = []
interpolated_y_values_1 = []
interpolated_y_values_2 = []


for x_val in x_values:
    y_val = evaluate_natural_cubic_spline(x_val, x_data, a, b, c, d)
    interpolated_y_values.append(y_val)

for x_val_1 in x_values_1:
    y_val_1 = evaluate_natural_cubic_spline(x_val_1, x_data_1, a1, b1, c1, d1)
    interpolated_y_values_1.append(y_val_1)

for x_val_2 in x_values_2:
    y_val_2 = evaluate_natural_cubic_spline(x_val_2, x_data_2, a2, b2, c2, d2)
    interpolated_y_values_2.append(y_val_2)


plt.figure(figsize=(10,10))
plt.plot(x_data, y_data, 'o', color='blue',  markersize=8)
plt.plot(x_data_1, y_data_1, 'o', color='red',  markersize=8)
plt.plot(x_data_2, y_data_2, 'o', color='green',  markersize=8)


plt.plot(x_values, interpolated_y_values, '-', color='blue', linewidth=2)
plt.plot(x_values_1, interpolated_y_values_1, '-', color='red', linewidth=2)
plt.plot(x_values_2, interpolated_y_values_2, '-', color='green', linewidth=2)

plt.title('Letter N')
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.grid(True)
plt.axis('equal')
plt.show()
