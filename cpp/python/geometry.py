# geometry.py
import numpy as np


class SolidAngle(object):
    def __init__(self, angle_range=[0, np.pi, -np.pi / 2, np.pi / 2], theta_inc=np.pi / 180, beta_inc=np.pi / 180,
                 unit="rad"):
        """

        :param angle_range: 
        :param theta_inc: 
        :param beta_inc: 
        :param unit: 
        """
        if unit == "deg":
            angle_range = np.radians(angle_range)
            theta_inc = np.radians(theta_inc)
            beta_inc = np.radians(beta_inc)
        self.range = angle_range
        self.theta_inc = theta_inc
        self.beta_inc = beta_inc
        self.theta_arr = np.arange(angle_range[0], angle_range[1], self.theta_inc)
        self.beta_arr = np.arange(angle_range[2], angle_range[3], self.beta_inc)
        self.subtend_calc = np.sum(np.sin(self.theta_arr) * self.theta_inc * self.beta_inc) * self.beta_arr.size
        self.subtend = (np.cos(self.theta_arr[0]) - np.cos(self.theta_arr[-1])) * (self.beta_arr[-1] - self.beta_arr[0])

    def shift(self, new_center):
        pass

    def slice(self, angle_range=[0, np.pi, -np.pi / 2, np.pi / 2]):
        """


        :rtype : 2D array.
        :param angle_range: 
        :return: 
        """
        theta_slice = np.logical_and(self.theta_arr >= angle_range[0], self.theta_arr <= angle_range[1])
        beta_slice = np.logical_and(self.beta_arr >= angle_range[2], self.beta_arr <= angle_range[3])
        omega_slice = np.multiply(*np.meshgrid(theta_slice, beta_slice))
        return omega_slice

    def domega(self, theta_rad):
        return (np.cos(theta_rad - self.theta_inc / 2.) - np.cos(theta_rad + self.theta_inc / 2.)) * self.beta_inc


class Illumination(object):
    def __init__(self, ev0=1e4, psi=0, alpha=0):
        self.ev0 = ev0
        self.psi = psi
        self.alpha = alpha
