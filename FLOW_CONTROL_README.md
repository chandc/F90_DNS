# DNS Flow Control Documentation

## Overview
This DNS solver now supports two different flow control methods to maintain specific flow conditions during simulation:

1. **Constant Pressure Gradient (Method 1)**: Maintains a fixed pressure gradient throughout the simulation
2. **Constant Volume Flow (Method 2)**: Uses a PI controller to dynamically adjust pressure gradient to maintain constant bulk velocity

## Configuration

Flow control parameters are specified in the input file under the `&flow_control` namelist:

```fortran
&flow_control
  flow_control_method = 1,           ! Method: 1=constant pressure, 2=constant volume
  target_pressure_gradient = 1.0,   ! Target pressure gradient (used in method 1)
  target_bulk_velocity = 1.0,       ! Target bulk velocity (used in method 2)
  controller_gain = 0.1,             ! PI controller gain (used in method 2)
  controller_update_freq = 10        ! Controller update frequency in time steps (used in method 2)
/
```

## Flow Control Methods

### Method 1: Constant Pressure Gradient
- **Use case**: When you want to maintain a fixed driving force
- **Implementation**: Applies a constant pressure gradient at each time step
- **Parameters**: 
  - `target_pressure_gradient`: The constant pressure gradient to apply
- **Advantages**: Simple, predictable, computationally efficient
- **Example**: `input_constant_pressure.dat`

### Method 2: Constant Volume Flow
- **Use case**: When you want to maintain a specific flow rate
- **Implementation**: PI controller adjusts pressure gradient to maintain target bulk velocity
- **Parameters**:
  - `target_bulk_velocity`: Desired bulk velocity to maintain
  - `controller_gain`: PI controller gain (typically 0.01-0.1)
  - `controller_update_freq`: How often to update controller (every N time steps)
- **Advantages**: Maintains precise flow rate, useful for flow rate-dependent studies
- **Example**: `input_constant_volume.dat`

## PI Controller Details

The PI controller for constant volume flow uses:
- **Proportional term**: Immediate response to velocity error
- **Integral term**: Eliminates steady-state error over time
- **Anti-windup**: Prevents controller saturation
- **Safety limits**: Pressure gradient bounded between 0.1 and 10.0

Controller equation:
```
dp/dt = Kp * error + Ki * integral_error
```

Where:
- `error = target_bulk_velocity - current_bulk_velocity`
- `Kp = Ki = controller_gain`

## Bulk Velocity Calculation

Bulk velocity is calculated using LGL quadrature integration:
```
U_bulk = (1/H) * ∫[0 to H] u(y) dy
```

Where:
- H is the channel height
- Integration uses LGL weights for spectral accuracy
- Calculation performed in physical space for accuracy

## Output and Monitoring

Flow control statistics are output during simulation including:
- Current pressure gradient
- Current bulk velocity (for method 2)
- Controller error and integral (for method 2)
- Update frequency information

## Usage Examples

### Constant Pressure Gradient
```bash
# Use input_constant_pressure.dat
./dns_solver input_constant_pressure.dat
```

### Constant Volume Flow
```bash
# Use input_constant_volume.dat  
./dns_solver input_constant_volume.dat
```

## Parameter Guidelines

### Controller Gain (`controller_gain`)
- Start with 0.05-0.1 for most cases
- Lower values (0.01-0.05): More stable, slower response
- Higher values (0.1-0.5): Faster response, potential instability
- Adjust based on flow Reynolds number and time step

### Update Frequency (`controller_update_freq`)
- Typical range: 5-20 time steps
- Lower values: More responsive control, higher computational cost
- Higher values: More stable, less responsive
- Should be much smaller than characteristic flow time scales

### Target Bulk Velocity (`target_bulk_velocity`)
- Should be physically reasonable for your Reynolds number
- For channel flow, typical values: 0.5-2.0 times laminar centerline velocity
- Monitor initial transients to ensure achievable target

## Troubleshooting

### Controller Oscillations
- Reduce `controller_gain`
- Increase `controller_update_freq`
- Check that target velocity is achievable

### Slow Controller Response
- Increase `controller_gain` (carefully)
- Reduce `controller_update_freq`
- Verify LGL grid resolution is adequate

### Pressure Gradient Limits
- If hitting safety limits (0.1 or 10.0), adjust target velocity
- Check for numerical instabilities
- Verify grid resolution and time step

## Implementation Notes

- Flow control is applied after computing the RHS but before the pressure projection
- LGL integration ensures spectral accuracy for bulk velocity calculation
- Controller state is maintained between time steps
- Default behavior (method 1 with gradient 1.0) maintains backward compatibility
