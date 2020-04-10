## Method 1 
 
- Wave equation: KDV 
$$
\dfrac{du}{dt} + 6u\dfrac{du}{dx}  + \dfrac{d^{3}u}{dx} = 0
$$

where u= u(x,.t) = velocity 


Simplify equation: 

$$
\dfrac{du}{dt} + 3\dfrac{du^2}{dx} + \dfrac{d^{3}u}{dx} = 0
$$

- Using **central difference formula** for **spatial** discretization:

$$
3\dfrac{du^2}{dx} = 3\big[ \dfrac{(u^2_{i+1})-(u^2_{i-1})}{2x}\big]
$$

$$
\dfrac{d^{3}u}{dx} = \big[ \dfrac{(u_{i+2})-(2u_{i+1})+(2u_{i-1})-(u_{i-2})}{2x^3}\big]
$$

- Using **Forward Euler Method ** for time discretization

$$
t^{N+1} = t^N + \Delta t
$$

$$
u_t = \dfrac{u_i^{N+1}-u_i^N}{\Delta t}
$$

Hence

$$
\dfrac{u_i^{N+1}-u_i^N}{\Delta t} + 3\big[ \dfrac{(u^2_{i+1})-(u^2_{i-1})}{2x}\big] + \big[ \dfrac{(u_{i+2})-(2u_{i+1})+(2u_{i-1})-(u_{i-2})}{2x^3}\big] = 0
$$


Final simplified equation: 
>$$
u_i^{N+1}= u_i^N-\Delta  t\cdot F(u^N)
$$

Actual solution: 
>$$
u(x,0 ) = \dfrac{1}{2}\sec{(h^2\dfrac{(x-t)}{2})} 
$$


space interval : [-10,10] , $\Delta x = 1$

time interval : [-10,10], $\Delta t = 0.01s$

## Method 2
- Wave equation: KDV  

- Using **central difference formula** for **spatial** discretization

- Using **third order Runge-Kutta Method** for time discretization:
$$
u^{N+1}= \dfrac{1}{3}u^N+\dfrac{2}{3}u^{(2)}+\dfrac{2}{3}\Delta t \cdot F(u^{(2)})
$$

$$
u^{(2)}= \dfrac{3}{4}u^N+\dfrac{1}{4}u^{(1)}+\dfrac{1}{4}\Delta t \cdot F(u^{(1)})
$$

$$
u^{(1)}= u^N+\Delta t \cdot F(u^{(N)})
$$



## Method 3
- Wave equation: Camassa-Holm (CH)
  $$
    m_t + um_x +2mu_x  = 0
  $$

 
  $$ m = u-u_{xx} $$

  where: $u(x,t)$ = velocity, 
$m(x,t)$ = momentum 

 Simplify equation: 
 $$
m_t +(um)_x +mu_x = 0 
$$
- Using **central difference formula** for **spatial** discretization
$$
(um)_x = \big[ \dfrac{(u_{i+1})-(u_{i-1}  \cdot m_{i-1})}{2x}\big]
$$

$$
mu_x = m_i\big[ \dfrac{(u_{i+1})-(u_{i-1})}{2x}\big]
$$

Rewritten equation:

>$$
m_t + \big[ \dfrac{(u_{i+1})-(u_{i-1}  \cdot m_{i-1})}{2x}\big] + m_i\big[ \dfrac{(u_{i+1})-(u_{i-1})}{2x}\big] = 0
$$

- Using **third order Runge-Kutta Method** for time discretization:
$$
u^{N+1}= \dfrac{1}{3}u^N+\dfrac{2}{3}u^{(2)}+\dfrac{2}{3}\Delta t \cdot F(u^{(2)})
$$

$$
u^{(2)}= \dfrac{3}{4}u^N+\dfrac{1}{4}u^{(1)}+\dfrac{1}{4}\Delta t \cdot F(u^{(1)})
$$

$$
u^{(1)}= u^N+\Delta t \cdot F(u^{(N)})
$$

- Discritize $m= u-u_{xx}$
$$
m_i = u_i - \big[ \dfrac{u_{i_2}-2u_i+u_{i-2}}{4x^3}\big]
$$

- solve linear algebra 




## Method 4 

- Wave equation: 2 component Camassa-Holm (CH)
  $$
    m_t + um_x +2mu_x  = - g\rho\rho_x
  $$

 
  $$ m = u-u_{xx} $$

  $$ \rho_t + (\rho u )_x = 0 $$

  where: 
  
  $u(x,t)$ = velocity, 
  
  $m(x,t)$ = momentum ,

  $\rho(x,t)$ = density, 

  g= gravity constant

 Simplified equation: 
 $$
m_t +(um)_x +mu_x = -\dfrac{g}{x}\rho^2x
$$
- Using **central difference formula** for **spatial** discretization
$$
(um)_x = \big[ \dfrac{(u_{i+1})-(u_{i-1}  \cdot m_{i-1})}{2x}\big]
$$

$$
mu_x = m_i\big[ \dfrac{(u_{i+1})-(u_{i-1})}{2x}\big]
$$

$$
-\dfrac{g}{x}\rho^2x = -\dfrac{g}{2}\big[ \dfrac{(\rho_{i+1}^2)-( \rho^2_{i-1})}{2x}\big]
$$

$$
(\rho u ) _x= \big[ \dfrac{(\rho_{i-1}* u_{i-1})-( \rho_{i+1}*u_{i+1})}{2x}\big]
$$
Rewritten equation:

>$$
m_t + \big[ \dfrac{(u_{i+1})-(u_{i-1}  \cdot m_{i-1})}{2x}\big] + m_i\big[ \dfrac{(u_{i+1})-(u_{i-1})}{2x}\big] = 0
$$

- Using **third order Runge-Kutta Method** for time discretization:
$$
u^{N+1}= \dfrac{1}{3}u^N+\dfrac{2}{3}u^{(2)}+\dfrac{2}{3}\Delta t \cdot F(u^{(2)})
$$

$$
u^{(2)}= \dfrac{3}{4}u^N+\dfrac{1}{4}u^{(1)}+\dfrac{1}{4}\Delta t \cdot F(u^{(1)})
$$

$$
u^{(1)}= u^N+\Delta t \cdot F(u^{(N)})
$$

- Discritize $m= u-u_{xx}$
$$
m_i = u_i - \big[ \dfrac{u_{i_2}-2u_i+u_{i-2}}{4x^3}\big]
$$

- solve linear algebra 
- initial condition (Dam Break Problem)
$$
\rho(x,0)=1+\tanh(x+4)-\tanh (x-4)  \\
U_0 = 0 
$$

- NO actual solution 


Add length scale-bifurication parameter ($\alpha$) to 
$$m=u-u_{xx} $$

So it becomes 

$$m=u-\alpha ^2u_{xx} $$

Try for $\alpha  = 0.1$

