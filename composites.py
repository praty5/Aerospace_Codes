import numpy as np
import matplotlib.pyplot as plt

# Material properties
E1 = 205  
E2 = 5   
E6 = 2.5 
v_1 = 0.25
v_2 = v_1*(E2/E1)

theta_degrees = np.arange(0, 91, 1)
theta_radians = np.deg2rad(theta_degrees)

Ex_values = []
Ey_values = []
Es_values = []
nu_xy_values = []
eta_xs_values = []
eta_ys_values = []

with open('orthotropic_properties.txt', 'w') as file:
    file.write("  Theta\t Ex (GPa)\t Ey (GPa)\t  Gxy(GPa)\t nu_xy\t  eta_xs\t  eta_ys\n")

    for theta in theta_radians:

        m = np.cos(theta)
        n = np.sin(theta)

        # Required quantities
        Ex = E1/(m**4 + (E1/E6 - 2*v_1)*(n**2)*(m**2) + E1/E2*(n**4))
         
        Ey = E2/(m**4 + (E2/E6 - 2*v_2)*(n**2)*(m**2) + E2/E1*(n**4))

        Es = E6/(n**4 + m**4 + 2*(2*(E6/E1)*(1+2*v_1)+2*E6/E2 -1)*(n**2)*(m**2))

        nu_xy = Ex*(m**2/E1*(m**2*v_1-n**2) + n**2/E2*(n**2*v_2-m**2) + m**2*n**2/E6)

        eta_xs = (2*m*n*(m**2 - (n**2)*v_1) - 2*m*n*E1/E2*(n**2 - (m**2)*v_2) + (E1/E6)*(m*n**3 - m**3*n)) / (m**4 + (E1/E6 - 2*v_1)*(n**2)*(m**2) + E1/E2*(n**4))
        
        eta_ys = ((-2)*m*n*(m**2 - (n**2)*v_2) + 2*m*n*E2/E1*(n**2 - (m**2)*v_1) + (E2/E6)*(m**3*n - m*n**3)) / (m**4 + (E2/E6 - 2*v_2)*(n**2)*(m**2) + E2/E1*(n**4))

        # Append the results to the lists
        Ex_values.append(Ex)
        Ey_values.append(Ey)
        Es_values.append(Es)
        nu_xy_values.append(nu_xy)
        eta_xs_values.append(eta_xs)
        eta_ys_values.append(eta_ys)

        file.write(f"{np.rad2deg(theta):>7.0f}\t{Ex:>8.4f}\t{Ey:>8.4f}\t{Es:>8.4f}\t{nu_xy:>8.4f}\t{eta_xs:>8.4f}\t{eta_ys:>8.4f}\n")

# Create plots
plt.figure(figsize=(15, 10))

plt.subplot(2, 3, 1)
plt.plot(theta_degrees, Ex_values)
plt.title('Ex vs. Theta')
plt.xlabel('Theta (degrees)')
plt.ylabel('Ex (GPa)')

plt.subplot(2, 3, 2)
plt.plot(theta_degrees, Ey_values)
plt.title('Ey vs. Theta')
plt.xlabel('Theta (degrees)')
plt.ylabel('Ey (GPa)')

plt.subplot(2, 3, 3)
plt.plot(theta_degrees, Es_values)
plt.title('Gxy vs. Theta')
plt.xlabel('Theta (degrees)')
plt.ylabel('Gxy (GPa)')

plt.subplot(2, 3, 4)
plt.plot(theta_degrees, nu_xy_values)
plt.title('nu_xy vs. Theta')
plt.xlabel('Theta (degrees)')
plt.ylabel('nu_xy')

plt.subplot(2, 3, 5)
plt.plot(theta_degrees, eta_xs_values)
plt.title('eta_xs vs. Theta')
plt.xlabel('Theta (degrees)')
plt.ylabel('eta_xs')

plt.subplot(2, 3, 6)
plt.plot(theta_degrees, eta_ys_values)
plt.title('eta_ys vs. Theta')
plt.xlabel('Theta (degrees)')
plt.ylabel('eta_ys')

plt.tight_layout()
plt.savefig('orthotropic_properties.png')

# Show the plot
plt.show()
