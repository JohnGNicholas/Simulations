#iterate over the range of values of interest
def irange(start, stop, step):
    while start < stop:
        yield round(start, 4)
        start += step

#Runge-Kutta method
def rk(odes, interval, step_size, initial_values):
    if len(odes) == len(initial_values):
        if len(odes) == 1:
            u = [initial_values]
            index = 0
            for t in interval[0:len(interval)-1]:
                u.append(u[index] + step_size * odes(t, u[index]))
                index += 1
        else:
            u = [[initial_values[0]]]
            for n in range(1, len(odes)):
                u[0].append(initial_values[n])
            index = 0
            for t in interval[0:len(interval)-1]:
                # Calculate y values for k1
                uValues1 = [u[index][0]]
                for n in range(1,len(odes)):
                    uValues1.append(u[index][n])
                # Calculate first slope 
                k1 = [odes[0](t, uValues1)]
                for n in range(1, len(odes)):
                    k1.append(odes[n](t, uValues1))
                # Calculate y values for k2
                uValues2 = [u[index][0] + 0.5 * step_size * k1[0]]
                for n in range(1, len(odes)):
                    uValues2.append(u[index][n] + 0.5 * step_size * k1[n])
                # Calculate second slope
                k2 = [odes[0](t + 0.5 * step_size, uValues2)]
                for n in range(1, len(odes)):
                    k2.append(odes[n](t + 0.5 * step_size, uValues2))
                # Calculate y values for k3
                uValues3 = [u[index][0] + 0.5 * step_size * k2[0]]
                for n in range(1, len(odes)):
                    uValues3.append(u[index][n] + 0.5 * step_size * k2[n])
                # Calculate third slope
                k3 = [odes[0](t + 0.5 * step_size, uValues3)]
                for n in range(1, len(odes)):
                    k3.append(odes[n](t + 0.5 * step_size, uValues3))
                # Calculate y values for fourth slope
                uValues4 = [u[index][0] + step_size * k3[0]]
                for n in range(1, len(odes)):
                    uValues4.append(u[index][n] + step_size * k3[n])
                # Calculate fourth slope
                k4 = [odes[0](t + step_size, uValues4)]
                for n in range(1, len(odes)):
                    k4.append(odes[n](t + step_size, uValues4))
                # Average slopes
                kbar = [(k1[0] + 2*k2[0] + 2*k3[0] + k4[0])/6]
                for n in range(1, len(odes)):
                    kbar.append((k1[n] + 2*k2[n] + 2*k3[n] + k4[n])/6)
                # Update Values in multidimensional list u
                uFinalValues = [u[index][0] + step_size * kbar[0]]
                for n in range(1,len(odes)):
                    uFinalValues.append(u[index][n] + step_size * kbar[n])
                u.append(uFinalValues)
                index += 1
    else:
        print('Incorrect Number of ODEs and Initial Values')
    return u