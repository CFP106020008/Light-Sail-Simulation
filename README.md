# Light Sail Simulation

###### tags: `GitHub`
View this article on [HackMD](https://hackmd.io/@wxH3br6dSwK7CMUJlrGf-w/H1JRQcSHd).

## Physics of light sail

In this simulation, we assume the simplest light sail configuration: a sigle flat panel which can control it's attitude.

![](https://i.imgur.com/652QQrF.png)

We define $\vec{r}$ to be the position, $\vec{v}$ to be the velocity and $\hat{p}$ to be the pointing of the sail.

Let's also define the acceleration of the sail which is facing the sun $(\theta=0)$ at $r=1AU$ as $a_E$.

We can see the acceleration of the sail is:

$$
\frac{d^2\vec{r}}{dt^2} = \vec{a} = \left(\frac{a_Er_E^2}{r^2}\cos{\theta}\right)\hat{p}$$

where

$$
\cos{\theta} = \frac{\vec{r}\cdot\hat{p}}{r}
$$

Then we can use ```scipy.integrate.solve_ivp``` to solve this initial value problem.


## Attitude Control
Since we want to accelerate our sail as fast as possible, we want the $\vec{a}\cdot\vec{v}$ to be as large as possible. The extrema should happens at:

$$
\frac{d(\vec{a}\cdot\vec{v})}{d\theta_p} = 0,~~ \theta_p = \tan^{-1}{(p_y/p_x)}
$$

Solving this equation we get:

$$
\hat{p} = (\frac{\vec{r}}{r} + \frac{\vec{v}}{v})/\sqrt{2},~~(\textrm{for acceleration})
$$

$$
\hat{p} = (\frac{\vec{r}}{r} - \frac{\vec{v}}{v})/\sqrt{2},~~(\textrm{for deceleration})
$$


Note that this is a native choose. In the real orbital design of light sails this should be more complicated. I control the pointing of the sail in the ```Decided_Pointing``` function.


## Initial Condition
In the code there are two choice of initial condition.

The first one is cicular, which put the sail in an circular orbit with configurable radius $R$.

The second one is elliptical, which puts the sail in an orbit with Aphelion of $R_{init}$ and perihelion of $R_\odot$.

## Example result
https://user-images.githubusercontent.com/48315222/202633569-fd440964-25c1-4207-ae08-62d72458b714.mp4


