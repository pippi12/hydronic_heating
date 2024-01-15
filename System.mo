package System
  package MyComponents
model HEX_HP
  import SI = Modelica.Units.SI;
  replaceable package Medium = Buildings.Media.Water "Medium in the pipe";
  inner Modelica.Fluid.System system annotation(
    Placement(visible = true, transformation(origin = {-128, 88}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  // parameters
  parameter Real K = 150 "1st delay K";
  parameter Modelica.Units.SI.Time T = 500 "1st delay T [s]";
  parameter Modelica.Units.SI.Time L = 100 "wasted time L [s]";
  parameter Modelica.Units.SI.Temperature InitialTemp = 273.15 + 30 "Initial Temperature [K]";
  // parts
  Modelica.Fluid.Interfaces.FluidPort_a port_a(redeclare package Medium = Medium) annotation(
    Placement(visible = true, transformation(origin = {-120, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-120, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Fluid.Interfaces.FluidPort_b port_b(redeclare package Medium = Medium) annotation(
    Placement(visible = true, transformation(origin = {120, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {80, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Fluid.Vessels.ClosedVolume volume(redeclare package Medium = Medium, T_start = InitialTemp, V = 0.00001, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, massDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, use_HeatTransfer = true, use_T_start = true, use_portsData = false, nPorts = 2) annotation(
    Placement(transformation(origin = {58, 10}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Interfaces.RealInput freq annotation(
    Placement(visible = true, transformation(origin = {-126, 46}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-160, 118}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow annotation(
    Placement(transformation(origin = {8, 46}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Continuous.FirstOrder firstOrder(T = T, initType = Modelica.Blocks.Types.Init.InitialOutput, k = K, y_start = 0) annotation(
    Placement(transformation(origin = {-24, 46}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Nonlinear.FixedDelay fixedDelay(delayTime = L) annotation(
    Placement(transformation(origin = {-66, 46}, extent = {{-10, -10}, {10, 10}})));
equation
  connect(port_a, volume.ports[1]) annotation(
    Line(points = {{-120, 0}, {58, 0}}));
  connect(port_b, volume.ports[2]) annotation(
    Line(points = {{120, 0}, {58, 0}}));
  connect(prescribedHeatFlow.port, volume.heatPort) annotation(
    Line(points = {{18, 46}, {32, 46}, {32, 10}, {48, 10}}, color = {191, 0, 0}));
  connect(firstOrder.y, prescribedHeatFlow.Q_flow) annotation(
    Line(points = {{-13, 46}, {-3, 46}}, color = {0, 0, 127}));
  connect(freq, fixedDelay.u) annotation(
    Line(points = {{-126, 46}, {-78, 46}}, color = {0, 0, 127}));
  connect(fixedDelay.y, firstOrder.u) annotation(
    Line(points = {{-54, 46}, {-36, 46}}, color = {0, 0, 127}));
  annotation(
    Icon(graphics = {Rectangle(origin = {-20, 50}, fillColor = {77, 77, 77}, fillPattern = FillPattern.Solid, extent = {{-100, 90}, {100, -90}}), Text(origin = {-23, -60}, lineColor = {0, 0, 255}, extent = {{-85, 24}, {85, -24}}, textString = "%name"), Rectangle(origin = {-70, 0}, lineColor = {0, 0, 255}, fillColor = {0, 0, 255}, fillPattern = FillPattern.Solid, extent = {{-50, 4}, {50, -4}}), Rectangle(origin = {30, 0}, lineColor = {255, 0, 0}, fillColor = {255, 0, 0}, fillPattern = FillPattern.Solid, extent = {{-50, 4}, {50, -4}}), Text(origin = {-160, 155}, extent = {{-26, 27}, {26, -27}}, textString = "f"), Line(origin = {-130, 118}, points = {{-10, 0}, {10, 0}, {10, 0}}, thickness = 0.5)}, coordinateSystem(extent = {{-200, 180}, {100, -80}})),
    Diagram(coordinateSystem(extent = {{-140, -140}, {140, 140}})));
end HEX_HP;

    model roomCycle_bak
      replaceable package Medium = Buildings.Media.Water;
      inner Modelica.Fluid.System system annotation(
        Placement(transformation(origin = {-154, 102}, extent = {{-10, -10}, {10, 10}})));
      // Parameters
      // Initialization
      parameter Modelica.Units.SI.Temperature medium_initT = 273.15 + 30 "Water Initial Temperature [K]";
      parameter Modelica.Units.SI.Temperature rad_initT = 273.15 + 40;
      //
      parameter Modelica.Units.SI.SpecificHeatCapacity pip_c = 386 "Specific heat of pipe wall material [J/kg.K]";
      parameter Modelica.Units.SI.Density pip_rho = 8960 "Density of pipe wall material [kg/m3]";
      parameter Modelica.Units.SI.Length pip_dh = 0.025 "Inner diameter of pipe [m]";
      parameter Modelica.Units.SI.Length pip_thickness = 0.005 "Thickness of pipe wall [m]";
      parameter Modelica.Units.SI.Length pip_dIns = 0.005 "Thickness of pipe insulation [m]";
      parameter Modelica.Units.SI.ThermalConductivity pip_kIns = 1 "Heat conductivity of pipe insulation [W/m.K]";
      //
      final parameter Modelica.Units.SI.MassFlowRate m_flow_nominal = 6 "Mass flow rate [kg/s]";
      final parameter Modelica.Units.SI.PressureDifference dp_nominal = 4500 "Design pressure drop";
      // Components
      Buildings.Fluid.Movers.FlowControlled_m_flow pump_m_flow(redeclare package Medium = Medium, m_flow_nominal = m_flow_nominal, use_inputFilter = false, massFlowRates = {0, 0.5, 1}*m_flow_nominal, inputType = Buildings.Fluid.Types.InputType.Continuous, energyDynamics = Modelica.Fluid.Types.Dynamics.SteadyState, nominalValuesDefineDefaultPressureCurve = true, allowFlowReversal = false, T_start = medium_initT) "Pump with m_flow input" annotation(
        Placement(transformation(origin = {-74, 38}, extent = {{-10, -10}, {10, 10}})));
      Buildings.Fluid.HeatExchangers.Radiators.RadiatorEN442_2 rad(redeclare package Medium = Medium, Q_flow_nominal = 5000, TAir_nominal = 293.15, T_a_nominal = 353.15, T_b_nominal = 333.15, T_start = medium_initT, allowFlowReversal = false, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, fraRad = 0, p_start = 101325) annotation(
        Placement(visible = true, transformation(origin = {24, 6}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      Modelica.Blocks.Sources.Ramp ramp1(duration = 1, height = 0.167, offset = 0, startTime = 0) annotation(
        Placement(transformation(origin = {-88, 61}, extent = {{-7, -7}, {7, 7}})));
      Buildings.Fluid.FixedResistances.PlugFlowPipe pipe1(redeclare package Medium = Medium, T_start_in = medium_initT, T_start_out = medium_initT, allowFlowReversal = false, cPip = pip_c, dIns = pip_dIns, dh = pip_dh, kIns = pip_kIns, length = 2, rhoPip = pip_rho, thickness = pip_thickness, m_flow_nominal = m_flow_nominal) annotation(
        Placement(transformation(origin = {-100, 38}, extent = {{-10, -10}, {10, 10}})));
      Buildings.Fluid.FixedResistances.PlugFlowPipe pipe2(redeclare package Medium = Medium, T_start_in = medium_initT, T_start_out = medium_initT, allowFlowReversal = false, cPip = pip_c, dIns = pip_dIns, dh = pip_dh, kIns = pip_kIns, length = 2, rhoPip = pip_rho, thickness = pip_thickness, m_flow_nominal = m_flow_nominal) annotation(
        Placement(transformation(origin = {2, 38}, extent = {{-10, -10}, {10, 10}})));
      Buildings.Fluid.FixedResistances.PlugFlowPipe pipe3(redeclare package Medium = Medium, T_start_in = medium_initT, T_start_out = medium_initT, allowFlowReversal = false, cPip = pip_c, dIns = pip_dIns, dh = pip_dh, kIns = pip_kIns, length = 2, rhoPip = pip_rho, thickness = pip_thickness, m_flow_nominal = m_flow_nominal) annotation(
        Placement(transformation(origin = {2, -32}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
      Buildings.Fluid.FixedResistances.PlugFlowPipe pipe4(redeclare package Medium = Medium, T_start_in = medium_initT, T_start_out = medium_initT, allowFlowReversal = false, cPip = pip_c, dIns = pip_dIns, dh = pip_dh, kIns = pip_kIns, length = 2, rhoPip = pip_rho, thickness = pip_thickness, m_flow_nominal = m_flow_nominal) annotation(
        Placement(transformation(origin = {-100, -32}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
      Buildings.Fluid.FixedResistances.PlugFlowPipe pipe5(redeclare package Medium = Medium, T_start_in = medium_initT, T_start_out = medium_initT, allowFlowReversal = true, cPip = pip_c, dIns = pip_dIns, dh = pip_dh, kIns = pip_kIns, length = 2, rhoPip = pip_rho, thickness = pip_thickness, m_flow_nominal = m_flow_nominal) annotation(
        Placement(transformation(origin = {-50, -6}, extent = {{10, -10}, {-10, 10}}, rotation = 90)));
      Buildings.Fluid.Actuators.Valves.TwoWayLinear val(CvData = Buildings.Fluid.Types.CvTypes.OpPoint, redeclare package Medium = Medium, dpValve_nominal(displayUnit = "kPa") = dp_nominal, m_flow_nominal = m_flow_nominal, use_inputFilter = false) annotation(
        Placement(transformation(origin = {-50, 24}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      Modelica.Blocks.Sources.Constant const(k = 0) annotation(
        Placement(transformation(origin = {-20, 1}, extent = {{-7, 7}, {7, -7}}, rotation = 90)));
      // Interfaces
      Modelica.Fluid.Interfaces.FluidPort_a port_a(redeclare package Medium = Medium) annotation(
        Placement(visible = true, transformation(origin = {-130, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Fluid.Interfaces.FluidPort_b port_b(redeclare package Medium = Medium) annotation(
        Placement(visible = true, transformation(origin = {-130, -32}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a port_Con annotation(
        Placement(visible = true, transformation(origin = {68, 22}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a port_Rad annotation(
        Placement(visible = true, transformation(origin = {68, -12}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Buildings.HeatTransfer.Sources.FixedTemperature outdoor(T = 283.15) annotation(
        Placement(transformation(origin = {-152, 74}, extent = {{-8, -8}, {8, 8}})));
    equation
      connect(rad.heatPortCon, port_Con) annotation(
        Line(points = {{32, 8}, {46, 8}, {46, 22}, {68, 22}}, color = {191, 0, 0}));
      connect(rad.heatPortRad, port_Rad) annotation(
        Line(points = {{32, 4}, {46, 4}, {46, -12}, {68, -12}}, color = {191, 0, 0}));
      connect(port_a, pipe1.port_a) annotation(
        Line(points = {{-130, 38}, {-110, 38}}));
      connect(pipe1.port_b, pump_m_flow.port_a) annotation(
        Line(points = {{-90, 38}, {-84, 38}}, color = {0, 127, 255}));
      connect(pump_m_flow.port_b, pipe2.port_a) annotation(
        Line(points = {{-64, 38}, {-8, 38}}, color = {0, 127, 255}));
      connect(pipe2.port_b, rad.port_a) annotation(
        Line(points = {{12, 38}, {24, 38}, {24, 16}}, color = {0, 127, 255}));
      connect(ramp1.y, pump_m_flow.m_flow_in) annotation(
        Line(points = {{-80, 61}, {-75.3, 61}, {-75.3, 50}}, color = {0, 0, 127}));
      connect(pipe4.port_b, port_b) annotation(
        Line(points = {{-110, -32}, {-130, -32}}, color = {0, 127, 255}));
      connect(pipe4.port_a, pipe3.port_b) annotation(
        Line(points = {{-90, -32}, {-8, -32}}, color = {0, 127, 255}));
      connect(pipe3.port_a, rad.port_b) annotation(
        Line(points = {{12, -32}, {24, -32}, {24, -4}}, color = {0, 127, 255}));
      connect(const.y, val.y) annotation(
        Line(points = {{-20, 8.7}, {-20, 23.7}, {-35, 23.7}}, color = {0, 0, 127}));
      connect(outdoor.port, pipe1.heatPort) annotation(
        Line(points = {{-144, 74}, {-100, 74}, {-100, 48}}, color = {191, 0, 0}));
      connect(outdoor.port, pipe2.heatPort) annotation(
        Line(points = {{-144, 74}, {2, 74}, {2, 48}}, color = {191, 0, 0}));
      connect(outdoor.port, pipe4.heatPort) annotation(
        Line(points = {{-144, 74}, {-140, 74}, {-140, -14}, {-100, -14}, {-100, -22}}, color = {191, 0, 0}));
      connect(outdoor.port, pipe5.heatPort) annotation(
        Line(points = {{-144, 74}, {-140, 74}, {-140, -6}, {-60, -6}}, color = {191, 0, 0}));
      connect(outdoor.port, pipe3.heatPort) annotation(
        Line(points = {{-144, 74}, {-140, 74}, {-140, -14}, {2, -14}, {2, -22}}, color = {191, 0, 0}));
      connect(val.port_b, pipe5.port_a) annotation(
        Line(points = {{-50, 14}, {-50, 4}}, color = {0, 127, 255}));
      connect(pipe5.port_b, pipe4.port_a) annotation(
        Line(points = {{-50, -16}, {-50, -32}, {-90, -32}}, color = {0, 127, 255}));
      connect(pump_m_flow.port_b, val.port_a) annotation(
        Line(points = {{-64, 38}, {-50, 38}, {-50, 34}}, color = {0, 127, 255}));
      annotation(
        Icon(graphics = {Rectangle(fillColor = {255, 255, 255}, lineThickness = 1, extent = {{-100, 100}, {100, -100}}), Text(origin = {2, -1}, extent = {{-72, 37}, {72, -37}}, textString = "%name", textStyle = {TextStyle.Bold}), Ellipse(origin = {0, -4}, extent = {{-82, 78}, {82, -78}})}, coordinateSystem(extent = {{-100, -100}, {100, 100}})),
        Diagram(coordinateSystem(extent = {{-180, 120}, {80, -60}})));
    end roomCycle_bak;

    model roomCycle01
      replaceable package Medium = Buildings.Media.Water;
      inner Modelica.Fluid.System system annotation(
        Placement(transformation(origin = {-154, 102}, extent = {{-10, -10}, {10, 10}})));
      // Parameters
      // Initialization
      parameter Modelica.Units.SI.Temperature medium_initT = 273.15 + 30 "Water Initial Temperature [K]";
      parameter Modelica.Units.SI.Temperature rad_initT = 273.15 + 40;
      //
      parameter Modelica.Units.SI.SpecificHeatCapacity pip_c = 386 "Specific heat of pipe wall material [J/kg.K]";
      parameter Modelica.Units.SI.Density pip_rho = 8960 "Density of pipe wall material [kg/m3]";
      parameter Modelica.Units.SI.Length pip_dh = 0.025 "Inner diameter of pipe [m]";
      parameter Modelica.Units.SI.Length pip_thickness = 0.005 "Thickness of pipe wall [m]";
      parameter Modelica.Units.SI.Length pip_dIns = 0.005 "Thickness of pipe insulation [m]";
      parameter Modelica.Units.SI.ThermalConductivity pip_kIns = 1 "Heat conductivity of pipe insulation [W/m.K]";
      //
      final parameter Modelica.Units.SI.MassFlowRate m_flow_nominal = 6 "Mass flow rate [kg/s]";
      final parameter Modelica.Units.SI.PressureDifference dp_nominal = 4500 "Design pressure drop";
      // Components
      Buildings.Fluid.Movers.FlowControlled_m_flow pump_m_flow(redeclare package Medium = Medium, m_flow_nominal = m_flow_nominal, use_inputFilter = false, massFlowRates = {0, 0.5, 1}*m_flow_nominal, inputType = Buildings.Fluid.Types.InputType.Continuous, energyDynamics = Modelica.Fluid.Types.Dynamics.SteadyState, nominalValuesDefineDefaultPressureCurve = true, allowFlowReversal = false, T_start = medium_initT) "Pump with m_flow input" annotation(
        Placement(transformation(origin = {-30, 38}, extent = {{-10, -10}, {10, 10}})));
      Buildings.Fluid.HeatExchangers.Radiators.RadiatorEN442_2 rad(redeclare package Medium = Medium, Q_flow_nominal = 5000, TAir_nominal = 293.15, T_a_nominal = 353.15, T_b_nominal = 333.15, T_start = medium_initT, allowFlowReversal = false, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, fraRad = 0, p_start = 101325) annotation(
        Placement(visible = true, transformation(origin = {24, 6}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      Modelica.Blocks.Sources.Ramp ramp1(duration = 1, height = 0.167, offset = 0, startTime = 0) annotation(
        Placement(transformation(origin = {-61, 53}, extent = {{-7, -7}, {7, 7}})));
      Buildings.Fluid.FixedResistances.PlugFlowPipe pipe1(redeclare package Medium = Medium, T_start_in = medium_initT, T_start_out = medium_initT, allowFlowReversal = false, cPip = pip_c, dIns = pip_dIns, dh = pip_dh, kIns = pip_kIns, length = 2, rhoPip = pip_rho, thickness = pip_thickness, m_flow_nominal = m_flow_nominal) annotation(
        Placement(transformation(origin = {-100, 38}, extent = {{-10, -10}, {10, 10}})));
      Buildings.Fluid.FixedResistances.PlugFlowPipe pipe2(redeclare package Medium = Medium, T_start_in = medium_initT, T_start_out = medium_initT, allowFlowReversal = false, cPip = pip_c, dIns = pip_dIns, dh = pip_dh, kIns = pip_kIns, length = 2, rhoPip = pip_rho, thickness = pip_thickness, m_flow_nominal = m_flow_nominal) annotation(
        Placement(transformation(origin = {2, 38}, extent = {{-10, -10}, {10, 10}})));
      Buildings.Fluid.FixedResistances.PlugFlowPipe pipe3(redeclare package Medium = Medium, T_start_in = medium_initT, T_start_out = medium_initT, allowFlowReversal = false, cPip = pip_c, dIns = pip_dIns, dh = pip_dh, kIns = pip_kIns, length = 2, rhoPip = pip_rho, thickness = pip_thickness, m_flow_nominal = m_flow_nominal) annotation(
        Placement(transformation(origin = {2, -32}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
      Buildings.Fluid.FixedResistances.PlugFlowPipe pipe4(redeclare package Medium = Medium, T_start_in = medium_initT, T_start_out = medium_initT, allowFlowReversal = false, cPip = pip_c, dIns = pip_dIns, dh = pip_dh, kIns = pip_kIns, length = 2, rhoPip = pip_rho, thickness = pip_thickness, m_flow_nominal = m_flow_nominal) annotation(
        Placement(transformation(origin = {-100, -32}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
      // Interfaces
      Modelica.Fluid.Interfaces.FluidPort_a port_a(redeclare package Medium = Medium) annotation(
        Placement(visible = true, transformation(origin = {-130, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Fluid.Interfaces.FluidPort_b port_b(redeclare package Medium = Medium) annotation(
        Placement(visible = true, transformation(origin = {-130, -32}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a port_Con annotation(
        Placement(visible = true, transformation(origin = {68, 22}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a port_Rad annotation(
        Placement(visible = true, transformation(origin = {68, -12}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Buildings.HeatTransfer.Sources.FixedTemperature outdoor(T = 283.15) annotation(
        Placement(transformation(origin = {-152, 74}, extent = {{-8, -8}, {8, 8}})));
    equation
      connect(rad.heatPortCon, port_Con) annotation(
        Line(points = {{32, 8}, {46, 8}, {46, 22}, {68, 22}}, color = {191, 0, 0}));
      connect(rad.heatPortRad, port_Rad) annotation(
        Line(points = {{32, 4}, {46, 4}, {46, -12}, {68, -12}}, color = {191, 0, 0}));
      connect(port_a, pipe1.port_a) annotation(
        Line(points = {{-130, 38}, {-110, 38}}));
      connect(pipe1.port_b, pump_m_flow.port_a) annotation(
        Line(points = {{-90, 38}, {-40, 38}}, color = {0, 127, 255}));
      connect(pump_m_flow.port_b, pipe2.port_a) annotation(
        Line(points = {{-20, 38}, {-8, 38}}, color = {0, 127, 255}));
      connect(pipe2.port_b, rad.port_a) annotation(
        Line(points = {{12, 38}, {24, 38}, {24, 16}}, color = {0, 127, 255}));
      connect(ramp1.y, pump_m_flow.m_flow_in) annotation(
        Line(points = {{-53, 53}, {-30, 53}, {-30, 50}}, color = {0, 0, 127}));
      connect(pipe4.port_b, port_b) annotation(
        Line(points = {{-110, -32}, {-130, -32}}, color = {0, 127, 255}));
      connect(pipe4.port_a, pipe3.port_b) annotation(
        Line(points = {{-90, -32}, {-8, -32}}, color = {0, 127, 255}));
      connect(pipe3.port_a, rad.port_b) annotation(
        Line(points = {{12, -32}, {24, -32}, {24, -4}}, color = {0, 127, 255}));
      connect(outdoor.port, pipe1.heatPort) annotation(
        Line(points = {{-144, 74}, {-100, 74}, {-100, 48}}, color = {191, 0, 0}));
      connect(outdoor.port, pipe2.heatPort) annotation(
        Line(points = {{-144, 74}, {2, 74}, {2, 48}}, color = {191, 0, 0}));
      connect(outdoor.port, pipe4.heatPort) annotation(
        Line(points = {{-144, 74}, {-140, 74}, {-140, -14}, {-100, -14}, {-100, -22}}, color = {191, 0, 0}));
      connect(outdoor.port, pipe3.heatPort) annotation(
        Line(points = {{-144, 74}, {-140, 74}, {-140, -14}, {2, -14}, {2, -22}}, color = {191, 0, 0}));
      annotation(
        Icon(graphics = {Rectangle(fillColor = {255, 255, 255}, lineThickness = 1, extent = {{-100, 100}, {100, -100}}), Text(origin = {2, -1}, extent = {{-72, 37}, {72, -37}}, textString = "%name", textStyle = {TextStyle.Bold}), Ellipse(origin = {0, -4}, extent = {{-82, 78}, {82, -78}})}, coordinateSystem(extent = {{-100, -100}, {100, 100}})),
        Diagram(coordinateSystem(extent = {{-180, 120}, {80, -60}})));
    end roomCycle01;

    model roomCycle
      replaceable package Medium = Buildings.Media.Water;
      replaceable package MediumAir = Buildings.Media.Air;
      inner Modelica.Fluid.System system annotation(
        Placement(transformation(origin = {-154, 102}, extent = {{-10, -10}, {10, 10}})));
      // Parameters
      // Initialization
      parameter Modelica.Units.SI.Temperature medium_initT = 273.15 + 30 "Water Initial Temperature [K]";
      parameter Modelica.Units.SI.Temperature mediumRoomAir_initT = 273.15 + 30 "RoomAir Initial Temperature [K]";
      parameter Modelica.Units.SI.Temperature rad_initT = 273.15 + 40;
      //
      parameter Modelica.Units.SI.SpecificHeatCapacity pip_c = 386 "Specific heat of pipe wall material [J/kg.K]";
      parameter Modelica.Units.SI.Density pip_rho = 8960 "Density of pipe wall material [kg/m3]";
      parameter Modelica.Units.SI.Length pip_dh = 0.025 "Inner diameter of pipe [m]";
      parameter Modelica.Units.SI.Length pip_thickness = 0.005 "Thickness of pipe wall [m]";
      parameter Modelica.Units.SI.Length pip_dIns = 0.005 "Thickness of pipe insulation [m]";
      parameter Modelica.Units.SI.ThermalConductivity pip_kIns = 1 "Heat conductivity of pipe insulation [W/m.K]";
      parameter Modelica.Units.SI.Temperature amb_T = 273.15 - 5 "room initial temperature [K]";
      parameter Real eps_HEX = 0.8 "Rossnay HEX coefficient [-]";
      //
      final parameter Modelica.Units.SI.MassFlowRate m_flow_nominal = 6 "Mass flow rate [kg/s]";
      final parameter Modelica.Units.SI.PressureDifference dp_nominal = 4500 "Design pressure drop";
      // Components
      Buildings.Fluid.Movers.FlowControlled_m_flow pump_m_flow(redeclare package Medium = Medium, m_flow_nominal = m_flow_nominal, use_inputFilter = false, massFlowRates = {0, 0.5, 1}*m_flow_nominal, inputType = Buildings.Fluid.Types.InputType.Continuous, energyDynamics = Modelica.Fluid.Types.Dynamics.SteadyState, nominalValuesDefineDefaultPressureCurve = true, allowFlowReversal = false, T_start = medium_initT) "Pump with m_flow input" annotation(
        Placement(transformation(origin = {-74, 38}, extent = {{-10, -10}, {10, 10}})));
      Buildings.Fluid.HeatExchangers.Radiators.RadiatorEN442_2 rad(redeclare package Medium = Medium, Q_flow_nominal = 5000, TAir_nominal = 293.15, T_a_nominal = 353.15, T_b_nominal = 333.15, T_start = medium_initT, allowFlowReversal = false, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, fraRad = 0, p_start = 101325) annotation(
        Placement(visible = true, transformation(origin = {24, 6}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      Modelica.Blocks.Sources.Ramp ramp1(duration = 1, height = 0.167, offset = 0, startTime = 0) annotation(
        Placement(transformation(origin = {-87, 62}, extent = {{-8, -8}, {8, 8}})));
      Buildings.Fluid.FixedResistances.PlugFlowPipe pipe1(redeclare package Medium = Medium, T_start_in = medium_initT, T_start_out = medium_initT, allowFlowReversal = false, cPip = pip_c, dIns = pip_dIns, dh = pip_dh, kIns = pip_kIns, length = 2, rhoPip = pip_rho, thickness = pip_thickness, m_flow_nominal = m_flow_nominal) annotation(
        Placement(transformation(origin = {-100, 38}, extent = {{-10, -10}, {10, 10}})));
      Buildings.Fluid.FixedResistances.PlugFlowPipe pipe2(redeclare package Medium = Medium, T_start_in = medium_initT, T_start_out = medium_initT, allowFlowReversal = false, cPip = pip_c, dIns = pip_dIns, dh = pip_dh, kIns = pip_kIns, length = 2, rhoPip = pip_rho, thickness = pip_thickness, m_flow_nominal = m_flow_nominal) annotation(
        Placement(transformation(origin = {2, 38}, extent = {{-10, -10}, {10, 10}})));
      Buildings.Fluid.FixedResistances.PlugFlowPipe pipe3(redeclare package Medium = Medium, T_start_in = medium_initT, T_start_out = medium_initT, allowFlowReversal = false, cPip = pip_c, dIns = pip_dIns, dh = pip_dh, kIns = pip_kIns, length = 2, rhoPip = pip_rho, thickness = pip_thickness, m_flow_nominal = m_flow_nominal) annotation(
        Placement(transformation(origin = {2, -32}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
      Buildings.Fluid.FixedResistances.PlugFlowPipe pipe4(redeclare package Medium = Medium, T_start_in = medium_initT, T_start_out = medium_initT, allowFlowReversal = false, cPip = pip_c, dIns = pip_dIns, dh = pip_dh, kIns = pip_kIns, length = 2, rhoPip = pip_rho, thickness = pip_thickness, m_flow_nominal = m_flow_nominal) annotation(
        Placement(transformation(origin = {-100, -32}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
      Buildings.Fluid.FixedResistances.PlugFlowPipe pipe5(redeclare package Medium = Medium, T_start_in = medium_initT, T_start_out = medium_initT, allowFlowReversal = true, cPip = pip_c, dIns = pip_dIns, dh = pip_dh, kIns = pip_kIns, length = 2, rhoPip = pip_rho, thickness = pip_thickness, m_flow_nominal = m_flow_nominal) annotation(
        Placement(transformation(origin = {-50, -6}, extent = {{10, -10}, {-10, 10}}, rotation = 90)));
      Buildings.Fluid.Actuators.Valves.TwoWayLinear val(CvData = Buildings.Fluid.Types.CvTypes.OpPoint, redeclare package Medium = Medium, dpValve_nominal(displayUnit = "kPa") = dp_nominal, m_flow_nominal = m_flow_nominal, use_inputFilter = false) annotation(
        Placement(transformation(origin = {-50, 24}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      Modelica.Blocks.Sources.Constant const(k = 0) annotation(
        Placement(transformation(origin = {-20, 1}, extent = {{-7, 7}, {7, -7}}, rotation = 90)));
      Buildings.Fluid.HeatExchangers.ConstantEffectiveness Lossnay(redeclare package Medium1 = MediumAir, redeclare package Medium2 = MediumAir, allowFlowReversal1 = false, allowFlowReversal2 = false, dp1_nominal = 500, dp2_nominal = 10, eps = eps_HEX, m1_flow_nominal = 5, m2_flow_nominal = 5) annotation(
        Placement(transformation(origin = {80, -32}, extent = {{-10, -10}, {10, 10}})));
      // Interfaces
      Modelica.Fluid.Interfaces.FluidPort_a port_a(redeclare package Medium = Medium) annotation(
        Placement(visible = true, transformation(origin = {-130, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Fluid.Interfaces.FluidPort_b port_b(redeclare package Medium = Medium) annotation(
        Placement(visible = true, transformation(origin = {-130, -32}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      //
      Buildings.HeatTransfer.Sources.FixedTemperature ambT(T = amb_T) annotation(
        Placement(transformation(origin = {-152, 74}, extent = {{-8, -8}, {8, 8}})));
      Buildings.Fluid.MixingVolumes.MixingVolume roomAir(redeclare package Medium = MediumAir, m_flow_nominal = 0.01, V = 10, p_start(displayUnit = "Pa") = 102315, T_start = 273.15 + 20, nPorts = 2, massDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial) annotation(
        Placement(transformation(origin = {58, 6}, extent = {{-10, -10}, {10, 10}})));
      Buildings.Fluid.FixedResistances.PlugFlowPipe pipeSA(redeclare package Medium = MediumAir, T_start_in = mediumRoomAir_initT, T_start_out = mediumRoomAir_initT, allowFlowReversal = true, cPip = pip_c, dIns = pip_dIns, dh = pip_dh, kIns = pip_kIns, length = 2, rhoPip = pip_rho, thickness = pip_thickness, m_flow_nominal = m_flow_nominal) annotation(
        Placement(transformation(origin = {120, -48}, extent = {{10, 10}, {-10, -10}})));
      Buildings.Fluid.FixedResistances.PlugFlowPipe pipeRA(redeclare package Medium = MediumAir, T_start_in = mediumRoomAir_initT, T_start_out = mediumRoomAir_initT, allowFlowReversal = false, cPip = pip_c, dIns = pip_dIns, dh = pip_dh, kIns = pip_kIns, length = 2, rhoPip = pip_rho, thickness = pip_thickness, m_flow_nominal = m_flow_nominal) annotation(
        Placement(transformation(origin = {174, -16}, extent = {{-10, -10}, {10, 10}})));
      Buildings.Fluid.Movers.FlowControlled_m_flow pump_m_flow_RA(redeclare package Medium = MediumAir, T_start = mediumRoomAir_initT, allowFlowReversal = true, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, inputType = Buildings.Fluid.Types.InputType.Continuous, m_flow_nominal = m_flow_nominal, massFlowRates = {0, 0.5, 1}*m_flow_nominal, nominalValuesDefineDefaultPressureCurve = true, use_inputFilter = false) annotation(
        Placement(transformation(origin = {132, -16}, extent = {{-10, -10}, {10, 10}})));
      Buildings.Fluid.Movers.FlowControlled_m_flow pump_m_flow_SA(redeclare package Medium = MediumAir, T_start = mediumRoomAir_initT, allowFlowReversal = false, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, inputType = Buildings.Fluid.Types.InputType.Continuous, m_flow_nominal = m_flow_nominal, massFlowRates = {0, 0.5, 1}*m_flow_nominal, nominalValuesDefineDefaultPressureCurve = true, use_inputFilter = false) annotation(
        Placement(transformation(origin = {164, -48}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
      Modelica.Blocks.Sources.Ramp ramp11(duration = 1, height = 1, offset = 0, startTime = 0) annotation(
        Placement(transformation(origin = {109, 12}, extent = {{-8, -8}, {8, 8}})));
      Buildings.Fluid.Sources.Boundary_pT ambAir(redeclare package Medium = MediumAir, p(displayUnit = "Pa") = 101325, T = 273.15 + 10, nPorts = 2) annotation(
        Placement(transformation(origin = {200, -16}, extent = {{10, -10}, {-10, 10}})));
    equation
      connect(port_a, pipe1.port_a) annotation(
        Line(points = {{-130, 38}, {-110, 38}}));
      connect(pipe1.port_b, pump_m_flow.port_a) annotation(
        Line(points = {{-90, 38}, {-84, 38}}, color = {0, 127, 255}));
      connect(pump_m_flow.port_b, pipe2.port_a) annotation(
        Line(points = {{-64, 38}, {-8, 38}}, color = {0, 127, 255}));
      connect(pipe2.port_b, rad.port_a) annotation(
        Line(points = {{12, 38}, {24, 38}, {24, 16}}, color = {0, 127, 255}));
      connect(ramp1.y, pump_m_flow.m_flow_in) annotation(
        Line(points = {{-78, 62}, {-75.3, 62}, {-75.3, 50}}, color = {0, 0, 127}));
      connect(pipe4.port_b, port_b) annotation(
        Line(points = {{-110, -32}, {-130, -32}}, color = {0, 127, 255}));
      connect(pipe4.port_a, pipe3.port_b) annotation(
        Line(points = {{-90, -32}, {-8, -32}}, color = {0, 127, 255}));
      connect(pipe3.port_a, rad.port_b) annotation(
        Line(points = {{12, -32}, {24, -32}, {24, -4}}, color = {0, 127, 255}));
      connect(const.y, val.y) annotation(
        Line(points = {{-20, 8.7}, {-20, 23.7}, {-35, 23.7}}, color = {0, 0, 127}));
      connect(val.port_b, pipe5.port_a) annotation(
        Line(points = {{-50, 14}, {-50, 4}}, color = {0, 127, 255}));
      connect(pipe5.port_b, pipe4.port_a) annotation(
        Line(points = {{-50, -16}, {-50, -32}, {-90, -32}}, color = {0, 127, 255}));
      connect(pump_m_flow.port_b, val.port_a) annotation(
        Line(points = {{-64, 38}, {-50, 38}, {-50, 34}}, color = {0, 127, 255}));
      connect(Lossnay.port_a2, pipeSA.port_b) annotation(
        Line(points = {{90, -38}, {105, -38}, {105, -48}, {110, -48}}, color = {0, 127, 255}));
      connect(pipeSA.port_a, pump_m_flow_SA.port_b) annotation(
        Line(points = {{130, -48}, {154, -48}}, color = {0, 127, 255}));
      connect(ramp11.y, pump_m_flow_RA.m_flow_in) annotation(
        Line(points = {{117.8, 12}, {131.8, 12}, {131.8, -4}}, color = {0, 0, 127}));
      connect(ramp11.y, pump_m_flow_SA.m_flow_in) annotation(
        Line(points = {{118, 12}, {162, 12}, {162, -36}, {164, -36}}, color = {0, 0, 127}));
      connect(roomAir.ports[1], Lossnay.port_a1) annotation(
        Line(points = {{58, -4}, {58, -26}, {70, -26}}, color = {0, 127, 255}));
      connect(roomAir.ports[2], Lossnay.port_b2) annotation(
        Line(points = {{58, -4}, {58, -38}, {70, -38}}, color = {0, 127, 255}));
      connect(rad.heatPortCon, roomAir.heatPort) annotation(
        Line(points = {{32, 8}, {48, 8}, {48, 6}}, color = {191, 0, 0}));
      connect(rad.heatPortRad, roomAir.heatPort) annotation(
        Line(points = {{32, 4}, {48, 4}, {48, 6}}, color = {191, 0, 0}));
      connect(pump_m_flow_SA.port_a, ambAir.ports[1]) annotation(
        Line(points = {{174, -48}, {190, -48}, {190, -16}}, color = {0, 127, 255}));
      connect(ambT.port, pipe1.heatPort) annotation(
        Line(points = {{-144, 74}, {-100, 74}, {-100, 48}}, color = {191, 0, 0}));
      connect(ambT.port, pipe2.heatPort) annotation(
        Line(points = {{-144, 74}, {2, 74}, {2, 48}}, color = {191, 0, 0}));
      connect(ambT.port, pipe5.heatPort) annotation(
        Line(points = {{-144, 74}, {-140, 74}, {-140, -6}, {-60, -6}}, color = {191, 0, 0}));
      connect(ambT.port, pipe4.heatPort) annotation(
        Line(points = {{-144, 74}, {-140, 74}, {-140, -6}, {-100, -6}, {-100, -22}}, color = {191, 0, 0}));
      connect(ambT.port, pipe3.heatPort) annotation(
        Line(points = {{-144, 74}, {-140, 74}, {-140, -6}, {-68, -6}, {-68, -22}, {2, -22}}, color = {191, 0, 0}));
      connect(ambT.port, pipeRA.heatPort) annotation(
        Line(points = {{-144, 74}, {120, 74}, {120, -6}, {174, -6}}, color = {191, 0, 0}));
      connect(ambT.port, pipeSA.heatPort) annotation(
        Line(points = {{-144, 74}, {-140, 74}, {-140, -58}, {120, -58}}, color = {191, 0, 0}));
      connect(pump_m_flow_RA.port_b, pipeRA.port_a) annotation(
        Line(points = {{142, -16}, {164, -16}}, color = {0, 127, 255}));
      connect(Lossnay.port_b1, pump_m_flow_RA.port_a) annotation(
        Line(points = {{90, -26}, {106, -26}, {106, -16}, {122, -16}}, color = {0, 127, 255}));
      connect(pipeRA.port_b, ambAir.ports[2]) annotation(
        Line(points = {{184, -16}, {190, -16}}, color = {0, 127, 255}));
      annotation(
        Icon(graphics = {Rectangle(fillColor = {255, 255, 255}, lineThickness = 1, extent = {{-100, 100}, {100, -100}}), Text(origin = {2, -1}, extent = {{-72, 37}, {72, -37}}, textString = "%name", textStyle = {TextStyle.Bold}), Ellipse(origin = {0, -4}, extent = {{-82, 78}, {82, -78}})}, coordinateSystem(extent = {{-100, -100}, {100, 100}})),
        Diagram(coordinateSystem(extent = {{-180, 120}, {220, -80}})));
    end roomCycle;

    model roomCycleWithoutLossnay_bak
      replaceable package Medium = Buildings.Media.Water;
      replaceable package MediumAir = Buildings.Media.Air;
      inner Modelica.Fluid.System system annotation(
        Placement(transformation(origin = {-154, 102}, extent = {{-10, -10}, {10, 10}})));
      // Parameters
      parameter Modelica.Units.SI.Temperature amb_T = 273.15 - 5 "Ambience initial temperature [K]";
      // Initialization
      parameter Modelica.Units.SI.Temperature medium_initT = 273.15 + 30 "Water Initial Temperature [K]";
      parameter Modelica.Units.SI.Temperature mediumRoomAir_initT = 273.15 + 30 "RoomAir Initial Temperature [K]";
      parameter Modelica.Units.SI.Pressure mediumRoomAir_initP = 101325 "RoomAir Initial Pressure [Pa]";
      parameter Modelica.Units.SI.Temperature rad_initT = 273.15 + 40;
      parameter Modelica.Units.SI.Temperature rw_initT = 273.15 + 30 "Room wall Initial Temperature [K]";
      // Pipe
      parameter Modelica.Units.SI.Length pip_len = 1 "Length of pipe wall [m]";
      parameter Modelica.Units.SI.Length pip5_len = 0.2 "Length of pipe wall [m]";
      parameter Modelica.Units.SI.SpecificHeatCapacity pip_c = 386 "Specific heat of pipe wall material [J/kg.K]";
      parameter Modelica.Units.SI.Density pip_rho = 8960 "Density of pipe wall material [kg/m3]";
      parameter Modelica.Units.SI.Length pip_dh = 0.025 "Inner diameter of pipe [m]";
      parameter Modelica.Units.SI.Length pip_thickness = 0.005 "Thickness of pipe wall [m]";
      parameter Modelica.Units.SI.Length pip_dIns = 0.005 "Thickness of pipe insulation [m]";
      parameter Modelica.Units.SI.ThermalConductivity pip_kIns = 398 "Heat conductivity of pipe insulation [W/m.K]";
      parameter Modelica.Units.SI.CoefficientOfHeatTransfer pip_air_htc = 10 "Air-Pipe heat transfer coeff. [W/(m2.K)]";
      // Room
      parameter Modelica.Units.SI.Length room_w = 1.6 "Width of room [m]";
      parameter Modelica.Units.SI.Length room_d = 4.8 "Depth of room [m]";
      parameter Modelica.Units.SI.Length room_h = 2.4 "Height of room [m]";
      parameter Modelica.Units.SI.Area room_area = 2*(room_w*room_d + room_d*room_h + room_h*room_w) "Area of room [m2]";
      parameter Modelica.Units.SI.Length rw_thickness = 0.1 "Thickness of room wall [m]";
      parameter Modelica.Units.SI.Density rw_rho = 114 "Density of room wall [m]";
      parameter Modelica.Units.SI.SpecificHeatCapacity rw_c = 1168 "Specific heat capacity of room wall [m]";
      parameter Modelica.Units.SI.ThermalConductivity rw_k = 0.3 "Heat conductivity of room wall [m]";
      parameter Modelica.Units.SI.Volume rw_vol = room_area*rw_thickness;
      parameter Modelica.Units.SI.CoefficientOfHeatTransfer rw_air_htc = 10 "Air-Room wall heat transfer coeff. [W/(m2.K)]";
      // Pump
      parameter Real m_flow = 5/60*0.001*1000 "mass flow [kg/s]";
      //
      final parameter Modelica.Units.SI.MassFlowRate m_flow_nominal = 6 "Mass flow rate [kg/s]";
      final parameter Modelica.Units.SI.PressureDifference dp_nominal = 4500 "Design pressure drop";
      // Components
      Buildings.Fluid.Movers.FlowControlled_m_flow pump_m_flow(redeclare package Medium = Medium, m_flow_nominal = m_flow_nominal, use_inputFilter = false, massFlowRates = {0, 0.5, 1}*m_flow_nominal, inputType = Buildings.Fluid.Types.InputType.Continuous, energyDynamics = Modelica.Fluid.Types.Dynamics.SteadyState, nominalValuesDefineDefaultPressureCurve = true, allowFlowReversal = false, T_start = medium_initT) "Pump with m_flow input" annotation(
        Placement(transformation(origin = {-74, 38}, extent = {{-10, -10}, {10, 10}})));
      Buildings.Fluid.HeatExchangers.Radiators.RadiatorEN442_2 rad(redeclare package Medium = Medium, Q_flow_nominal = 5000, TAir_nominal = 293.15, T_a_nominal = 353.15, T_b_nominal = 333.15, T_start = medium_initT, allowFlowReversal = false, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, fraRad = 0, p_start = 101325) annotation(
        Placement(visible = true, transformation(origin = {24, 6}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      Modelica.Blocks.Sources.Ramp ramp1(duration = 1, height = m_flow, offset = 0, startTime = 0) annotation(
        Placement(transformation(origin = {-87, 62}, extent = {{-8, -8}, {8, 8}})));
      Buildings.Fluid.FixedResistances.PlugFlowPipe pipe1(redeclare package Medium = Medium, T_start_in = medium_initT, T_start_out = medium_initT, allowFlowReversal = false, cPip = pip_c, dIns = pip_thickness, dh = pip_dh, kIns = pip_kIns, length = pip_len, rhoPip = pip_rho, thickness = pip_thickness, m_flow_nominal = m_flow_nominal, R = 1/(Modelica.Math.log((pip_dh + pip_thickness)/pip_dh)/(2*Modelica.Constants.pi*pip_kIns*pip_len) + 1/pip_air_htc/Modelica.Constants.pi/(pip_dh + pip_thickness)/pip_len)) annotation(
        Placement(transformation(origin = {-100, 38}, extent = {{-10, -10}, {10, 10}})));
      Buildings.Fluid.FixedResistances.PlugFlowPipe pipe2(redeclare package Medium = Medium, T_start_in = medium_initT, T_start_out = medium_initT, allowFlowReversal = false, cPip = pip_c, dIns = pip_thickness, dh = pip_dh, kIns = pip_kIns, length = pip_len, rhoPip = pip_rho, thickness = pip_thickness, m_flow_nominal = m_flow_nominal, R = 1/(Modelica.Math.log((pip_dh + pip_thickness)/pip_dh)/(2*Modelica.Constants.pi*pip_kIns*pip_len) + 1/pip_air_htc/Modelica.Constants.pi/(pip_dh + pip_thickness)/pip_len)) annotation(
        Placement(transformation(origin = {2, 38}, extent = {{-10, -10}, {10, 10}})));
      Buildings.Fluid.FixedResistances.PlugFlowPipe pipe3(redeclare package Medium = Medium, T_start_in = medium_initT, T_start_out = medium_initT, allowFlowReversal = false, cPip = pip_c, dIns = pip_thickness, dh = pip_dh, kIns = pip_kIns, length = pip_len, rhoPip = pip_rho, thickness = pip_thickness, m_flow_nominal = m_flow_nominal, R = 1/(Modelica.Math.log((pip_dh + pip_thickness)/pip_dh)/(2*Modelica.Constants.pi*pip_kIns*pip_len) + 1/pip_air_htc/Modelica.Constants.pi/(pip_dh + pip_thickness)/pip_len)) annotation(
        Placement(transformation(origin = {2, -32}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
      Buildings.Fluid.FixedResistances.PlugFlowPipe pipe4(redeclare package Medium = Medium, T_start_in = medium_initT, T_start_out = medium_initT, allowFlowReversal = false, cPip = pip_c, dIns = pip_thickness, dh = pip_dh, kIns = pip_kIns, length = pip_len, rhoPip = pip_rho, thickness = pip_thickness, m_flow_nominal = m_flow_nominal, R = 1/(Modelica.Math.log((pip_dh + pip_thickness)/pip_dh)/(2*Modelica.Constants.pi*pip_kIns*pip_len) + 1/pip_air_htc/Modelica.Constants.pi/(pip_dh + pip_thickness)/pip_len)) annotation(
        Placement(transformation(origin = {-100, -32}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
      Buildings.Fluid.FixedResistances.PlugFlowPipe pipe5(redeclare package Medium = Medium, T_start_in = medium_initT, T_start_out = medium_initT, allowFlowReversal = true, cPip = pip_c, dIns = pip_thickness, dh = pip_dh, kIns = pip_kIns, length = pip5_len, rhoPip = pip_rho, thickness = pip_thickness, m_flow_nominal = m_flow_nominal, R = 1/(Modelica.Math.log((pip_dh + pip_thickness)/pip_dh)/(2*Modelica.Constants.pi*pip_kIns*pip5_len) + 1/pip_air_htc/Modelica.Constants.pi/(pip_dh + pip_thickness)/pip5_len)) annotation(
        Placement(transformation(origin = {-50, -6}, extent = {{10, -10}, {-10, 10}}, rotation = 90)));
      Buildings.Fluid.Actuators.Valves.TwoWayLinear val(CvData = Buildings.Fluid.Types.CvTypes.OpPoint, redeclare package Medium = Medium, dpValve_nominal(displayUnit = "kPa") = dp_nominal, m_flow_nominal = m_flow_nominal, use_inputFilter = false) annotation(
        Placement(transformation(origin = {-50, 24}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      Modelica.Blocks.Sources.Constant const(k = 0) annotation(
        Placement(transformation(origin = {-20, 1}, extent = {{-7, 7}, {7, -7}}, rotation = 90)));
      // Interfaces
      Modelica.Fluid.Interfaces.FluidPort_a port_a(redeclare package Medium = Medium) annotation(
        Placement(visible = true, transformation(origin = {-130, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Fluid.Interfaces.FluidPort_b port_b(redeclare package Medium = Medium) annotation(
        Placement(visible = true, transformation(origin = {-130, -32}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      //
      Buildings.HeatTransfer.Sources.FixedTemperature ambT(T = amb_T) annotation(
        Placement(transformation(origin = {-152, 74}, extent = {{-8, -8}, {8, 8}})));
      Buildings.Fluid.MixingVolumes.MixingVolume roomAir(redeclare package Medium = MediumAir, m_flow_nominal = 0.01, V = room_w*room_d*room_h, p_start(displayUnit = "Pa") = mediumRoomAir_initP, T_start = mediumRoomAir_initT, nPorts = 0, massDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial) annotation(
        Placement(transformation(origin = {58, 6}, extent = {{-10, -10}, {10, 10}})));
      Modelica.Thermal.HeatTransfer.Components.HeatCapacitor heatCapacitor(C = rw_vol*rw_rho*rw_c, T(fixed = true, start = rw_initT)) annotation(
        Placement(transformation(origin = {106, 60}, extent = {{-10, -10}, {10, 10}})));
      Modelica.Thermal.HeatTransfer.Components.Convection conv_amb_wall annotation(
        Placement(transformation(origin = {48, 40}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
      Modelica.Blocks.Sources.Constant const_room(k = rw_air_htc*room_area) annotation(
        Placement(transformation(origin = {20, 90}, extent = {{-10, -10}, {10, 10}})));
      Modelica.Thermal.HeatTransfer.Components.Convection conv_room_wall annotation(
        Placement(transformation(origin = {164, 40}, extent = {{-10, -10}, {10, 10}})));
      Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor_wall1(G = rw_k*room_area*2/rw_thickness) annotation(
        Placement(transformation(origin = {82, 40}, extent = {{-10, -10}, {10, 10}})));
      Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor_wall2(G = rw_k*room_area*2/rw_thickness) annotation(
        Placement(transformation(origin = {130, 40}, extent = {{10, -10}, {-10, 10}})));
      Buildings.HeatTransfer.Convection.Exterior con(A = room_area, conMod = Buildings.HeatTransfer.Types.ExteriorConvection.Fixed, hFixed = rw_air_htc, azi = 1.570796326794897) annotation(
        Placement(transformation(origin = {86, -48}, extent = {{-10, -10}, {10, 10}})));
    equation
      connect(port_a, pipe1.port_a) annotation(
        Line(points = {{-130, 38}, {-110, 38}}));
      connect(pipe1.port_b, pump_m_flow.port_a) annotation(
        Line(points = {{-90, 38}, {-84, 38}}, color = {0, 127, 255}));
      connect(pump_m_flow.port_b, pipe2.port_a) annotation(
        Line(points = {{-64, 38}, {-8, 38}}, color = {0, 127, 255}));
      connect(pipe2.port_b, rad.port_a) annotation(
        Line(points = {{12, 38}, {24, 38}, {24, 16}}, color = {0, 127, 255}));
      connect(ramp1.y, pump_m_flow.m_flow_in) annotation(
        Line(points = {{-78, 62}, {-75.3, 62}, {-75.3, 50}}, color = {0, 0, 127}));
      connect(pipe4.port_b, port_b) annotation(
        Line(points = {{-110, -32}, {-130, -32}}, color = {0, 127, 255}));
      connect(pipe4.port_a, pipe3.port_b) annotation(
        Line(points = {{-90, -32}, {-8, -32}}, color = {0, 127, 255}));
      connect(pipe3.port_a, rad.port_b) annotation(
        Line(points = {{12, -32}, {24, -32}, {24, -4}}, color = {0, 127, 255}));
      connect(const.y, val.y) annotation(
        Line(points = {{-20, 8.7}, {-20, 23.7}, {-35, 23.7}}, color = {0, 0, 127}));
      connect(val.port_b, pipe5.port_a) annotation(
        Line(points = {{-50, 14}, {-50, 4}}, color = {0, 127, 255}));
      connect(pipe5.port_b, pipe4.port_a) annotation(
        Line(points = {{-50, -16}, {-50, -32}, {-90, -32}}, color = {0, 127, 255}));
      connect(pump_m_flow.port_b, val.port_a) annotation(
        Line(points = {{-64, 38}, {-50, 38}, {-50, 34}}, color = {0, 127, 255}));
      connect(rad.heatPortCon, roomAir.heatPort) annotation(
        Line(points = {{32, 8}, {48, 8}, {48, 6}}, color = {191, 0, 0}));
      connect(rad.heatPortRad, roomAir.heatPort) annotation(
        Line(points = {{32, 4}, {48, 4}, {48, 6}}, color = {191, 0, 0}));
      connect(ambT.port, pipe1.heatPort) annotation(
        Line(points = {{-144, 74}, {-100, 74}, {-100, 48}}, color = {191, 0, 0}));
      connect(ambT.port, pipe2.heatPort) annotation(
        Line(points = {{-144, 74}, {2, 74}, {2, 48}}, color = {191, 0, 0}));
      connect(ambT.port, pipe5.heatPort) annotation(
        Line(points = {{-144, 74}, {-140, 74}, {-140, -6}, {-60, -6}}, color = {191, 0, 0}));
      connect(ambT.port, pipe4.heatPort) annotation(
        Line(points = {{-144, 74}, {-140, 74}, {-140, -6}, {-100, -6}, {-100, -22}}, color = {191, 0, 0}));
      connect(ambT.port, pipe3.heatPort) annotation(
        Line(points = {{-144, 74}, {-140, 74}, {-140, -6}, {-68, -6}, {-68, -22}, {2, -22}}, color = {191, 0, 0}));
      connect(ambT.port, conv_amb_wall.fluid) annotation(
        Line(points = {{-144, 74}, {38, 74}, {38, 40}}, color = {191, 0, 0}));
      connect(conv_amb_wall.solid, thermalConductor_wall1.port_a) annotation(
        Line(points = {{58, 40}, {72, 40}}, color = {191, 0, 0}));
      connect(thermalConductor_wall1.port_b, heatCapacitor.port) annotation(
        Line(points = {{92, 40}, {100, 40}, {100, 50}, {106, 50}}, color = {191, 0, 0}));
      connect(heatCapacitor.port, thermalConductor_wall2.port_b) annotation(
        Line(points = {{106, 50}, {112, 50}, {112, 40}, {120, 40}}, color = {191, 0, 0}));
      connect(thermalConductor_wall2.port_a, conv_room_wall.solid) annotation(
        Line(points = {{140, 40}, {154, 40}}, color = {191, 0, 0}));
      connect(const_room.y, conv_amb_wall.Gc) annotation(
        Line(points = {{32, 90}, {48, 90}, {48, 50}}, color = {0, 0, 127}));
      connect(const_room.y, conv_room_wall.Gc) annotation(
        Line(points = {{32, 90}, {164, 90}, {164, 50}}, color = {0, 0, 127}));
      connect(roomAir.heatPort, conv_room_wall.fluid) annotation(
        Line(points = {{48, 6}, {40, 6}, {40, -10}, {182, -10}, {182, 40}, {174, 40}}, color = {191, 0, 0}));
      annotation(
        Icon(graphics = {Rectangle(fillColor = {255, 255, 255}, lineThickness = 1, extent = {{-100, 100}, {100, -100}}), Text(origin = {2, -1}, extent = {{-72, 37}, {72, -37}}, textString = "%name", textStyle = {TextStyle.Bold}), Ellipse(origin = {0, -4}, extent = {{-82, 78}, {82, -78}})}, coordinateSystem(extent = {{-100, -100}, {100, 100}})),
        Diagram(coordinateSystem(extent = {{-180, 120}, {220, -80}})));
    end roomCycleWithoutLossnay_bak;

    model roomCycleWithoutLossnay
      replaceable package Medium = Buildings.Media.Water;
      replaceable package MediumAir = Buildings.Media.Air;
      inner Modelica.Fluid.System system annotation(
        Placement(transformation(origin = {-154, 102}, extent = {{-10, -10}, {10, 10}})));
      // Parameters
      parameter Modelica.Units.SI.Temperature amb_T = 273.15 + 5 "Ambience initial temperature";
      // Initialization
      parameter Modelica.Units.SI.Temperature medium_initT = 273.15 + 30 "Water initial temperature";
      parameter Modelica.Units.SI.Temperature mediumRoomAir_initT = 273.15 + 20 "Room air initial temperature";
      parameter Modelica.Units.SI.Pressure mediumRoomAir_initP = 101325 "Room air initial Pressure";
      parameter Modelica.Units.SI.Temperature rad_initT = 273.15 + 40 "Radiator water initial temperature";
      parameter Modelica.Units.SI.Temperature rw_initT = 273.15 + 15 "Room wall initial temperature";
      // Pipe
      parameter Modelica.Units.SI.Length pip_len = 1 "Length of pipe wall";
      parameter Modelica.Units.SI.Length pip5_len = 0.2 "Length of pipe wall";
      parameter Modelica.Units.SI.SpecificHeatCapacity pip_c = 386 "Specific heat of pipe wall material";
      parameter Modelica.Units.SI.Density pip_rho = 8960 "Density of pipe wall material";
      parameter Modelica.Units.SI.Length pip_dh = 0.025 "Inner diameter of pipe";
      parameter Modelica.Units.SI.Length pip_thickness = 0.005 "Thickness of pipe wall";
      parameter Modelica.Units.SI.Length pip_dIns = 0.005 "Thickness of pipe insulation";
      parameter Modelica.Units.SI.ThermalConductivity pip_kIns = 398 "Heat conductivity of pipe insulation";
      parameter Modelica.Units.SI.CoefficientOfHeatTransfer pip_air_htc = 10 "Air-Pipe heat transfer coeff.";
      // Radiator
      parameter Modelica.Units.SI.Power q_flow_nominal = 5000 "Rated heat dissipation amount";
      // Room
      parameter Modelica.Units.SI.Length room_w = 1.6 "Width of room";
      parameter Modelica.Units.SI.Length room_d = 4.8 "Depth of room";
      parameter Modelica.Units.SI.Length room_h = 2.4 "Height of room";
      parameter Modelica.Units.SI.Area room_area = 2*(room_w*room_d + room_d*room_h + room_h*room_w) "Area of room";
      parameter Modelica.Units.SI.Length rw_thickness = 0.1 "Thickness of room wall";
      parameter Modelica.Units.SI.Density rw_rho = 144 "Density of room wall";
      parameter Modelica.Units.SI.SpecificHeatCapacity rw_c = 1168 "Specific heat capacity of room wall";
      parameter Modelica.Units.SI.ThermalConductivity rw_k = 0.3 "Heat conductivity of room wall";
      parameter Modelica.Units.SI.Volume rw_vol = room_area*rw_thickness;
      parameter Modelica.Units.SI.CoefficientOfHeatTransfer rw_air_htc = 5 "Air-Room wall heat transfer coeff.";
      // Pump
      parameter Real m_flow = 5/60*0.001*1000 "mass flow";
      //
      final parameter Modelica.Units.SI.MassFlowRate m_flow_nominal = 6 "Mass flow rate";
      final parameter Modelica.Units.SI.PressureDifference dp_nominal = 4500 "Design pressure drop";
      // Components
      Buildings.Fluid.Movers.FlowControlled_m_flow pump_m_flow(redeclare package Medium = Medium, m_flow_nominal = m_flow_nominal, use_inputFilter = false, massFlowRates = {0, 0.5, 1}*m_flow_nominal, inputType = Buildings.Fluid.Types.InputType.Continuous, energyDynamics = Modelica.Fluid.Types.Dynamics.SteadyState, nominalValuesDefineDefaultPressureCurve = true, allowFlowReversal = false, T_start = medium_initT) "Pump with m_flow input" annotation(
        Placement(transformation(origin = {-74, 38}, extent = {{-10, -10}, {10, 10}})));
      Buildings.Fluid.HeatExchangers.Radiators.RadiatorEN442_2 rad(redeclare package Medium = Medium, Q_flow_nominal = q_flow_nominal, TAir_nominal = 293.15, T_a_nominal = 353.15, T_b_nominal = 333.15, T_start = rad_initT, allowFlowReversal = false, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, fraRad = 0, p_start = 101325, nEle = 1) annotation(
        Placement(visible = true, transformation(origin = {24, 6}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      Modelica.Blocks.Sources.Ramp ramp1(duration = 1, height = m_flow, offset = 0, startTime = 0) annotation(
        Placement(transformation(origin = {-87, 62}, extent = {{-8, -8}, {8, 8}})));
      Buildings.Fluid.FixedResistances.PlugFlowPipe pipe1(redeclare package Medium = Medium, T_start_in = medium_initT, T_start_out = medium_initT, allowFlowReversal = false, cPip = pip_c, dIns = pip_thickness, dh = pip_dh, kIns = pip_kIns, length = pip_len, rhoPip = pip_rho, thickness = pip_thickness, m_flow_nominal = m_flow_nominal, R = Modelica.Math.log((pip_dh + pip_thickness)/pip_dh)/(2*Modelica.Constants.pi*pip_kIns) + 1/pip_air_htc/Modelica.Constants.pi/(pip_dh + pip_thickness)) annotation(
        Placement(transformation(origin = {-100, 38}, extent = {{-10, -10}, {10, 10}})));
      Buildings.Fluid.FixedResistances.PlugFlowPipe pipe2(redeclare package Medium = Medium, T_start_in = medium_initT, T_start_out = medium_initT, allowFlowReversal = false, cPip = pip_c, dIns = pip_thickness, dh = pip_dh, kIns = pip_kIns, length = pip_len, rhoPip = pip_rho, thickness = pip_thickness, m_flow_nominal = m_flow_nominal, R = Modelica.Math.log((pip_dh + pip_thickness)/pip_dh)/(2*Modelica.Constants.pi*pip_kIns) + 1/pip_air_htc/Modelica.Constants.pi/(pip_dh + pip_thickness)) annotation(
        Placement(transformation(origin = {2, 38}, extent = {{-10, -10}, {10, 10}})));
      Buildings.Fluid.FixedResistances.PlugFlowPipe pipe3(redeclare package Medium = Medium, T_start_in = medium_initT, T_start_out = medium_initT, allowFlowReversal = false, cPip = pip_c, dIns = pip_thickness, dh = pip_dh, kIns = pip_kIns, length = pip_len, rhoPip = pip_rho, thickness = pip_thickness, m_flow_nominal = m_flow_nominal, R = Modelica.Math.log((pip_dh + pip_thickness)/pip_dh)/(2*Modelica.Constants.pi*pip_kIns) + 1/pip_air_htc/Modelica.Constants.pi/(pip_dh + pip_thickness)) annotation(
        Placement(transformation(origin = {2, -32}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
      Buildings.Fluid.FixedResistances.PlugFlowPipe pipe4(redeclare package Medium = Medium, T_start_in = medium_initT, T_start_out = medium_initT, allowFlowReversal = false, cPip = pip_c, dIns = pip_thickness, dh = pip_dh, kIns = pip_kIns, length = pip_len, rhoPip = pip_rho, thickness = pip_thickness, m_flow_nominal = m_flow_nominal, R = Modelica.Math.log((pip_dh + pip_thickness)/pip_dh)/(2*Modelica.Constants.pi*pip_kIns) + 1/pip_air_htc/Modelica.Constants.pi/(pip_dh + pip_thickness)) annotation(
        Placement(transformation(origin = {-100, -32}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
      Buildings.Fluid.FixedResistances.PlugFlowPipe pipe5(redeclare package Medium = Medium, T_start_in = medium_initT, T_start_out = medium_initT, allowFlowReversal = true, cPip = pip_c, dIns = pip_thickness, dh = pip_dh, kIns = pip_kIns, length = pip5_len, rhoPip = pip_rho, thickness = pip_thickness, m_flow_nominal = m_flow_nominal, R = Modelica.Math.log((pip_dh + pip_thickness)/pip_dh)/(2*Modelica.Constants.pi*pip_kIns) + 1/pip_air_htc/Modelica.Constants.pi/(pip_dh + pip_thickness)) annotation(
        Placement(transformation(origin = {-50, -6}, extent = {{10, -10}, {-10, 10}}, rotation = 90)));
      Buildings.Fluid.Actuators.Valves.TwoWayLinear val(CvData = Buildings.Fluid.Types.CvTypes.OpPoint, redeclare package Medium = Medium, dpValve_nominal(displayUnit = "kPa") = dp_nominal, m_flow_nominal = m_flow_nominal, use_inputFilter = false) annotation(
        Placement(transformation(origin = {-50, 24}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      Modelica.Blocks.Sources.Constant const(k = 0) annotation(
        Placement(transformation(origin = {-20, 1}, extent = {{-7, 7}, {7, -7}}, rotation = 90)));
      // Interfaces
      Modelica.Fluid.Interfaces.FluidPort_a port_a(redeclare package Medium = Medium) annotation(
        Placement(visible = true, transformation(origin = {-130, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Fluid.Interfaces.FluidPort_b port_b(redeclare package Medium = Medium) annotation(
        Placement(visible = true, transformation(origin = {-130, -32}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      //
      Buildings.HeatTransfer.Sources.FixedTemperature ambT(T = amb_T) annotation(
        Placement(transformation(origin = {-152, 74}, extent = {{-8, -8}, {8, 8}})));
      Buildings.Fluid.MixingVolumes.MixingVolume roomAir(redeclare package Medium = MediumAir, m_flow_nominal = 0.01, V = room_w*room_d*room_h, p_start(displayUnit = "Pa") = mediumRoomAir_initP, T_start = mediumRoomAir_initT, nPorts = 0, massDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial) annotation(
        Placement(transformation(origin = {58, 6}, extent = {{-10, -10}, {10, 10}})));
      Modelica.Thermal.HeatTransfer.Components.HeatCapacitor heatCapacitor(C = rw_vol*rw_rho*rw_c, T(fixed = true, start = rw_initT)) annotation(
        Placement(transformation(origin = {106, 60}, extent = {{-10, -10}, {10, 10}})));
      Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor_wall1(G = rw_k*room_area*2/rw_thickness) annotation(
        Placement(transformation(origin = {82, 40}, extent = {{-10, -10}, {10, 10}})));
      Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor_wall2(G = rw_k*room_area*2/rw_thickness) annotation(
        Placement(transformation(origin = {130, 40}, extent = {{10, -10}, {-10, 10}})));
      Buildings.HeatTransfer.Convection.Exterior con_ext(A = room_area, conMod = Buildings.HeatTransfer.Types.ExteriorConvection.Fixed, hFixed = rw_air_htc, azi = 3.14/2) annotation(
        Placement(transformation(origin = {42, 40}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
      Buildings.HeatTransfer.Convection.Interior con_int(conMod = Buildings.HeatTransfer.Types.InteriorConvection.Fixed, A = room_area, hFixed = rw_air_htc) annotation(
        Placement(transformation(origin = {168, 40}, extent = {{-10, -10}, {10, 10}})));
      Modelica.Blocks.Sources.Constant const_v(k = 0) annotation(
        Placement(transformation(origin = {80, 90}, extent = {{-7, 7}, {7, -7}}, rotation = -180)));
      Modelica.Blocks.Sources.Constant const_d(k = 0) annotation(
        Placement(transformation(origin = {80, 70}, extent = {{-7, 7}, {7, -7}}, rotation = -180)));
    equation
      connect(port_a, pipe1.port_a) annotation(
        Line(points = {{-130, 38}, {-110, 38}}));
      connect(pipe1.port_b, pump_m_flow.port_a) annotation(
        Line(points = {{-90, 38}, {-84, 38}}, color = {0, 127, 255}));
      connect(pump_m_flow.port_b, pipe2.port_a) annotation(
        Line(points = {{-64, 38}, {-8, 38}}, color = {0, 127, 255}));
      connect(pipe2.port_b, rad.port_a) annotation(
        Line(points = {{12, 38}, {24, 38}, {24, 16}}, color = {0, 127, 255}));
      connect(ramp1.y, pump_m_flow.m_flow_in) annotation(
        Line(points = {{-78, 62}, {-75.3, 62}, {-75.3, 50}}, color = {0, 0, 127}));
      connect(pipe4.port_b, port_b) annotation(
        Line(points = {{-110, -32}, {-130, -32}}, color = {0, 127, 255}));
      connect(pipe4.port_a, pipe3.port_b) annotation(
        Line(points = {{-90, -32}, {-8, -32}}, color = {0, 127, 255}));
      connect(pipe3.port_a, rad.port_b) annotation(
        Line(points = {{12, -32}, {24, -32}, {24, -4}}, color = {0, 127, 255}));
      connect(const.y, val.y) annotation(
        Line(points = {{-20, 8.7}, {-20, 23.7}, {-35, 23.7}}, color = {0, 0, 127}));
      connect(val.port_b, pipe5.port_a) annotation(
        Line(points = {{-50, 14}, {-50, 4}}, color = {0, 127, 255}));
      connect(pipe5.port_b, pipe4.port_a) annotation(
        Line(points = {{-50, -16}, {-50, -32}, {-90, -32}}, color = {0, 127, 255}));
      connect(pump_m_flow.port_b, val.port_a) annotation(
        Line(points = {{-64, 38}, {-50, 38}, {-50, 34}}, color = {0, 127, 255}));
      connect(rad.heatPortCon, roomAir.heatPort) annotation(
        Line(points = {{32, 8}, {48, 8}, {48, 6}}, color = {191, 0, 0}));
      connect(rad.heatPortRad, roomAir.heatPort) annotation(
        Line(points = {{32, 4}, {48, 4}, {48, 6}}, color = {191, 0, 0}));
      connect(ambT.port, pipe1.heatPort) annotation(
        Line(points = {{-144, 74}, {-100, 74}, {-100, 48}}, color = {191, 0, 0}));
      connect(ambT.port, pipe2.heatPort) annotation(
        Line(points = {{-144, 74}, {2, 74}, {2, 48}}, color = {191, 0, 0}));
      connect(ambT.port, pipe5.heatPort) annotation(
        Line(points = {{-144, 74}, {-140, 74}, {-140, -6}, {-60, -6}}, color = {191, 0, 0}));
      connect(ambT.port, pipe4.heatPort) annotation(
        Line(points = {{-144, 74}, {-140, 74}, {-140, -6}, {-100, -6}, {-100, -22}}, color = {191, 0, 0}));
      connect(ambT.port, pipe3.heatPort) annotation(
        Line(points = {{-144, 74}, {-140, 74}, {-140, -6}, {-68, -6}, {-68, -22}, {2, -22}}, color = {191, 0, 0}));
      connect(thermalConductor_wall1.port_b, heatCapacitor.port) annotation(
        Line(points = {{92, 40}, {100, 40}, {100, 50}, {106, 50}}, color = {191, 0, 0}));
      connect(heatCapacitor.port, thermalConductor_wall2.port_b) annotation(
        Line(points = {{106, 50}, {112, 50}, {112, 40}, {120, 40}}, color = {191, 0, 0}));
      connect(con_ext.solid, thermalConductor_wall1.port_a) annotation(
        Line(points = {{52, 40}, {72, 40}}, color = {191, 0, 0}));
      connect(con_ext.fluid, ambT.port) annotation(
        Line(points = {{32, 40}, {20, 40}, {20, 74}, {-144, 74}}, color = {191, 0, 0}));
      connect(thermalConductor_wall2.port_a, con_int.solid) annotation(
        Line(points = {{140, 40}, {158, 40}}, color = {191, 0, 0}));
      connect(con_int.fluid, roomAir.heatPort) annotation(
        Line(points = {{178, 40}, {182, 40}, {182, -10}, {48, -10}, {48, 6}}, color = {191, 0, 0}));
      connect(const_v.y, con_ext.v) annotation(
        Line(points = {{72, 90}, {62, 90}, {62, 50}, {54, 50}}, color = {0, 0, 127}));
      connect(const_d.y, con_ext.dir) annotation(
        Line(points = {{72, 70}, {64, 70}, {64, 46}, {54, 46}}, color = {0, 0, 127}));
      annotation(
        Icon(graphics = {Rectangle(fillColor = {255, 255, 255}, lineThickness = 1, extent = {{-100, 100}, {100, -100}}), Text(origin = {2, -1}, extent = {{-72, 37}, {72, -37}}, textString = "%name", textStyle = {TextStyle.Bold}), Ellipse(origin = {0, -4}, extent = {{-82, 78}, {82, -78}})}, coordinateSystem(extent = {{-100, -100}, {100, 100}})),
        Diagram(coordinateSystem(extent = {{-180, 120}, {220, -80}})));
    end roomCycleWithoutLossnay;

    model Lossnay
      replaceable package MediumAir = Buildings.Media.Air;
      inner Modelica.Fluid.System system annotation(
        Placement(transformation(origin = {12, 36}, extent = {{-10, -10}, {10, 10}})));
      // Parameters
      // Initialization
      parameter Modelica.Units.SI.Temperature mediumRoomAir_initT = 273.15 + 30 "RoomAir Initial Temperature [K]";
      //
      parameter Modelica.Units.SI.SpecificHeatCapacity pip_c = 386 "Specific heat of pipe wall material [J/kg.K]";
      parameter Modelica.Units.SI.Density pip_rho = 8960 "Density of pipe wall material [kg/m3]";
      parameter Modelica.Units.SI.Length pip_dh = 0.025 "Inner diameter of pipe [m]";
      parameter Modelica.Units.SI.Length pip_thickness = 0.005 "Thickness of pipe wall [m]";
      parameter Modelica.Units.SI.Length pip_dIns = 0.005 "Thickness of pipe insulation [m]";
      parameter Modelica.Units.SI.ThermalConductivity pip_kIns = 1 "Heat conductivity of pipe insulation [W/m.K]";
      parameter Modelica.Units.SI.Temperature amb_T = 273.15 - 5 "room initial temperature [K]";
      parameter Real eps_HEX = 0.8 "Rossnay HEX coefficient [-]";
      //
      final parameter Modelica.Units.SI.MassFlowRate m_flow_nominal = 6 "Mass flow rate [kg/s]";
      final parameter Modelica.Units.SI.PressureDifference dp_nominal = 4500 "Design pressure drop";
      // Components
      Buildings.Fluid.HeatExchangers.ConstantEffectiveness Lossnay(redeclare package Medium1 = MediumAir, redeclare package Medium2 = MediumAir, allowFlowReversal1 = false, allowFlowReversal2 = false, dp1_nominal = 10, dp2_nominal = 10, eps = eps_HEX, m1_flow_nominal = 5, m2_flow_nominal = 5) annotation(
        Placement(transformation(origin = {80, -32}, extent = {{-10, -10}, {10, 10}})));
      // Interfaces
      Modelica.Fluid.Interfaces.FluidPort_a port_a(redeclare package Medium = MediumAir) annotation(
        Placement(transformation(origin = {48, -18}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {-100, 50}, extent = {{-10, -10}, {10, 10}})));
      Modelica.Fluid.Interfaces.FluidPort_b port_b(redeclare package Medium = MediumAir) annotation(
        Placement(transformation(origin = {48, -48}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {-100, -50}, extent = {{-10, -10}, {10, 10}})));
      //
      Buildings.HeatTransfer.Sources.FixedTemperature ambT(T = amb_T) annotation(
        Placement(transformation(origin = {10, 8}, extent = {{-8, -8}, {8, 8}})));
      Buildings.Fluid.FixedResistances.PlugFlowPipe pipeSA(redeclare package Medium = MediumAir, T_start_in = mediumRoomAir_initT, T_start_out = mediumRoomAir_initT, allowFlowReversal = true, cPip = pip_c, dIns = pip_dIns, dh = pip_dh, kIns = pip_kIns, length = 2, rhoPip = pip_rho, thickness = pip_thickness, m_flow_nominal = m_flow_nominal) annotation(
        Placement(transformation(origin = {120, -48}, extent = {{10, 10}, {-10, -10}})));
      Buildings.Fluid.FixedResistances.PlugFlowPipe pipeRA(redeclare package Medium = MediumAir, T_start_in = mediumRoomAir_initT, T_start_out = mediumRoomAir_initT, allowFlowReversal = false, cPip = pip_c, dIns = pip_dIns, dh = pip_dh, kIns = pip_kIns, length = 2, rhoPip = pip_rho, thickness = pip_thickness, m_flow_nominal = m_flow_nominal) annotation(
        Placement(transformation(origin = {174, -16}, extent = {{-10, -10}, {10, 10}})));
      Buildings.Fluid.Movers.FlowControlled_m_flow pump_m_flow_SA(redeclare package Medium = MediumAir, T_start = mediumRoomAir_initT, allowFlowReversal = true, energyDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial, inputType = Buildings.Fluid.Types.InputType.Continuous, m_flow_nominal = m_flow_nominal, massFlowRates = {0, 0.5, 1}*m_flow_nominal, nominalValuesDefineDefaultPressureCurve = true, use_inputFilter = false) annotation(
        Placement(transformation(origin = {164, -48}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
      Modelica.Blocks.Sources.Ramp ramp11(duration = 1, height = 0.5, offset = 0.5, startTime = 0) annotation(
        Placement(transformation(origin = {109, 30}, extent = {{-8, -8}, {8, 8}})));
      Buildings.Fluid.Sources.Boundary_pT ambAir(redeclare package Medium = MediumAir, p(displayUnit = "Pa") = 101325, T = 273.15 + 10, nPorts = 2) annotation(
        Placement(transformation(origin = {200, -16}, extent = {{10, -10}, {-10, 10}})));
    equation
      connect(Lossnay.port_a2, pipeSA.port_b) annotation(
        Line(points = {{90, -38}, {105, -38}, {105, -48}, {110, -48}}, color = {0, 127, 255}));
      connect(pipeSA.port_a, pump_m_flow_SA.port_b) annotation(
        Line(points = {{130, -48}, {154, -48}}, color = {0, 127, 255}));
      connect(ramp11.y, pump_m_flow_SA.m_flow_in) annotation(
        Line(points = {{118, 30}, {162, 30}, {162, -36}, {164, -36}}, color = {0, 0, 127}));
      connect(pump_m_flow_SA.port_a, ambAir.ports[1]) annotation(
        Line(points = {{174, -48}, {190, -48}, {190, -16}}, color = {0, 127, 255}));
      connect(ambT.port, pipeRA.heatPort) annotation(
        Line(points = {{18, 8}, {120, 8}, {120, -6}, {174, -6}}, color = {191, 0, 0}));
      connect(ambT.port, pipeSA.heatPort) annotation(
        Line(points = {{18, 8}, {30, 8}, {30, -58}, {120, -58}}, color = {191, 0, 0}));
      connect(pipeRA.port_b, ambAir.ports[2]) annotation(
        Line(points = {{184, -16}, {190, -16}}, color = {0, 127, 255}));
      connect(Lossnay.port_a1, port_a) annotation(
        Line(points = {{70, -26}, {60, -26}, {60, -18}, {48, -18}}, color = {0, 127, 255}));
      connect(Lossnay.port_b2, port_b) annotation(
        Line(points = {{70, -38}, {60, -38}, {60, -48}, {48, -48}}, color = {0, 127, 255}));
      connect(Lossnay.port_b1, pipeRA.port_a) annotation(
        Line(points = {{90, -26}, {140, -26}, {140, -16}, {164, -16}}, color = {0, 127, 255}));
      annotation(
        Icon(graphics = {Rectangle(fillColor = {255, 255, 255}, lineThickness = 1, extent = {{-100, 100}, {100, -100}}), Text(origin = {-10, 120}, extent = {{-72, 37}, {72, -37}}, textString = "%name", textStyle = {TextStyle.Bold}), Ellipse(origin = {0, -4}, extent = {{-82, 78}, {82, -78}})}, coordinateSystem(extent = {{-100, -100}, {100, 100}})),
        Diagram(coordinateSystem(extent = {{0, 60}, {220, -60}})));
    end Lossnay;

    model Lossnay_bak
      replaceable package MediumAir = Buildings.Media.Air;
      inner Modelica.Fluid.System system annotation(
        Placement(transformation(origin = {12, 36}, extent = {{-10, -10}, {10, 10}})));
      // Parameters
      // Initialization
      parameter Modelica.Units.SI.Temperature mediumRoomAir_initT = 273.15 + 30 "RoomAir Initial Temperature [K]";
      //
      parameter Modelica.Units.SI.SpecificHeatCapacity pip_c = 386 "Specific heat of pipe wall material [J/kg.K]";
      parameter Modelica.Units.SI.Density pip_rho = 8960 "Density of pipe wall material [kg/m3]";
      parameter Modelica.Units.SI.Length pip_dh = 0.025 "Inner diameter of pipe [m]";
      parameter Modelica.Units.SI.Length pip_thickness = 0.005 "Thickness of pipe wall [m]";
      parameter Modelica.Units.SI.Length pip_dIns = 0.005 "Thickness of pipe insulation [m]";
      parameter Modelica.Units.SI.ThermalConductivity pip_kIns = 1 "Heat conductivity of pipe insulation [W/m.K]";
      parameter Modelica.Units.SI.Temperature amb_T = 273.15 - 5 "room initial temperature [K]";
      parameter Real eps_HEX = 0.8 "Rossnay HEX coefficient [-]";
      //
      final parameter Modelica.Units.SI.MassFlowRate m_flow_nominal = 6 "Mass flow rate [kg/s]";
      final parameter Modelica.Units.SI.PressureDifference dp_nominal = 4500 "Design pressure drop";
      // Components
      Buildings.Fluid.HeatExchangers.ConstantEffectiveness Lossnay(redeclare package Medium1 = MediumAir, redeclare package Medium2 = MediumAir, allowFlowReversal1 = false, allowFlowReversal2 = false, dp1_nominal = 10, dp2_nominal = 10, eps = eps_HEX, m1_flow_nominal = 5, m2_flow_nominal = 5) annotation(
        Placement(transformation(origin = {80, -32}, extent = {{-10, -10}, {10, 10}})));
      // Interfaces
      Modelica.Fluid.Interfaces.FluidPort_a port_a(redeclare package Medium = MediumAir) annotation(
        Placement(transformation(origin = {48, -18}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {-100, 50}, extent = {{-10, -10}, {10, 10}})));
      Modelica.Fluid.Interfaces.FluidPort_b port_b(redeclare package Medium = MediumAir) annotation(
        Placement(transformation(origin = {48, -48}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {-100, -50}, extent = {{-10, -10}, {10, 10}})));
      //
      Buildings.HeatTransfer.Sources.FixedTemperature ambT(T = amb_T) annotation(
        Placement(transformation(origin = {10, 8}, extent = {{-8, -8}, {8, 8}})));
      Buildings.Fluid.FixedResistances.PlugFlowPipe pipeSA(redeclare package Medium = MediumAir, T_start_in = mediumRoomAir_initT, T_start_out = mediumRoomAir_initT, allowFlowReversal = true, cPip = pip_c, dIns = pip_dIns, dh = pip_dh, kIns = pip_kIns, length = 2, rhoPip = pip_rho, thickness = pip_thickness, m_flow_nominal = m_flow_nominal) annotation(
        Placement(transformation(origin = {120, -48}, extent = {{10, 10}, {-10, -10}})));
      Buildings.Fluid.FixedResistances.PlugFlowPipe pipeRA(redeclare package Medium = MediumAir, T_start_in = mediumRoomAir_initT, T_start_out = mediumRoomAir_initT, allowFlowReversal = false, cPip = pip_c, dIns = pip_dIns, dh = pip_dh, kIns = pip_kIns, length = 2, rhoPip = pip_rho, thickness = pip_thickness, m_flow_nominal = m_flow_nominal) annotation(
        Placement(transformation(origin = {174, -16}, extent = {{-10, -10}, {10, 10}})));
      Buildings.Fluid.Movers.FlowControlled_m_flow pump_m_flow_RA(redeclare package Medium = MediumAir, T_start = mediumRoomAir_initT, allowFlowReversal = true, energyDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial, inputType = Buildings.Fluid.Types.InputType.Continuous, m_flow_nominal = m_flow_nominal, massFlowRates = {0, 0.5, 1}*m_flow_nominal, nominalValuesDefineDefaultPressureCurve = true, use_inputFilter = false) annotation(
        Placement(transformation(origin = {132, -16}, extent = {{-10, -10}, {10, 10}})));
      Buildings.Fluid.Movers.FlowControlled_m_flow pump_m_flow_SA(redeclare package Medium = MediumAir, T_start = mediumRoomAir_initT, allowFlowReversal = true, energyDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial, inputType = Buildings.Fluid.Types.InputType.Continuous, m_flow_nominal = m_flow_nominal, massFlowRates = {0, 0.5, 1}*m_flow_nominal, nominalValuesDefineDefaultPressureCurve = true, use_inputFilter = false) annotation(
        Placement(transformation(origin = {164, -48}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
      Modelica.Blocks.Sources.Ramp ramp11(duration = 1, height = 0.5, offset = 0.5, startTime = 0) annotation(
        Placement(transformation(origin = {109, 30}, extent = {{-8, -8}, {8, 8}})));
      Buildings.Fluid.Sources.Boundary_pT ambAir(redeclare package Medium = MediumAir, p(displayUnit = "Pa") = 101325, T = 273.15 + 10, nPorts = 2) annotation(
        Placement(transformation(origin = {200, -16}, extent = {{10, -10}, {-10, 10}})));
    equation
      connect(Lossnay.port_a2, pipeSA.port_b) annotation(
        Line(points = {{90, -38}, {105, -38}, {105, -48}, {110, -48}}, color = {0, 127, 255}));
      connect(pipeSA.port_a, pump_m_flow_SA.port_b) annotation(
        Line(points = {{130, -48}, {154, -48}}, color = {0, 127, 255}));
      connect(ramp11.y, pump_m_flow_RA.m_flow_in) annotation(
        Line(points = {{118, 30}, {118, -4}, {131.8, -4}}, color = {0, 0, 127}));
      connect(ramp11.y, pump_m_flow_SA.m_flow_in) annotation(
        Line(points = {{118, 30}, {162, 30}, {162, -36}, {164, -36}}, color = {0, 0, 127}));
      connect(pump_m_flow_SA.port_a, ambAir.ports[1]) annotation(
        Line(points = {{174, -48}, {190, -48}, {190, -16}}, color = {0, 127, 255}));
      connect(ambT.port, pipeRA.heatPort) annotation(
        Line(points = {{18, 8}, {120, 8}, {120, -6}, {174, -6}}, color = {191, 0, 0}));
      connect(ambT.port, pipeSA.heatPort) annotation(
        Line(points = {{18, 8}, {30, 8}, {30, -58}, {120, -58}}, color = {191, 0, 0}));
      connect(pump_m_flow_RA.port_b, pipeRA.port_a) annotation(
        Line(points = {{142, -16}, {164, -16}}, color = {0, 127, 255}));
      connect(Lossnay.port_b1, pump_m_flow_RA.port_a) annotation(
        Line(points = {{90, -26}, {106, -26}, {106, -16}, {122, -16}}, color = {0, 127, 255}));
      connect(pipeRA.port_b, ambAir.ports[2]) annotation(
        Line(points = {{184, -16}, {190, -16}}, color = {0, 127, 255}));
      connect(Lossnay.port_a1, port_a) annotation(
        Line(points = {{70, -26}, {60, -26}, {60, -18}, {48, -18}}, color = {0, 127, 255}));
      connect(Lossnay.port_b2, port_b) annotation(
        Line(points = {{70, -38}, {60, -38}, {60, -48}, {48, -48}}, color = {0, 127, 255}));
      annotation(
        Icon(graphics = {Rectangle(fillColor = {255, 255, 255}, lineThickness = 1, extent = {{-100, 100}, {100, -100}}), Text(origin = {-10, 120}, extent = {{-72, 37}, {72, -37}}, textString = "%name", textStyle = {TextStyle.Bold}), Ellipse(origin = {0, -4}, extent = {{-82, 78}, {82, -78}})}, coordinateSystem(extent = {{-100, -100}, {100, 100}})),
        Diagram(coordinateSystem(extent = {{0, 60}, {220, -60}})));
    end Lossnay_bak;

    model HexControl
      //parameters
      parameter Modelica.Units.SI.Real dT = 1.5 "Temperature difference from OFF to ON";
      parameter Modelica.Units.SI.Real lower = 20 "Lower limit frequency";
      parameter Modelica.Units.SI.Real upper = 80 "Upper limit frequency";
    
    Modelica.Blocks.Interfaces.RealInput T_set annotation(
        Placement(visible = true, transformation(origin = {-102, 42}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-100, 60}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
    Modelica.Blocks.Interfaces.RealInput T_BT annotation(
        Placement(visible = true, transformation(origin = {-104, -28}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-100, -60}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
    Modelica.Blocks.Nonlinear.Limiter limiterMax(uMax = upper, uMin = 0)  annotation(
        Placement(visible = true, transformation(origin = {34, 42}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Math.Feedback feedback1 annotation(
        Placement(visible = true, transformation(origin = {-40, 42}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Continuous.PI pi(T = 500, initType = Modelica.Blocks.Types.Init.InitialState, k = 5, x_start = 0, y_start = 0)  annotation(
        Placement(visible = true, transformation(origin = {-2, 42}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Logical.Switch switch1 annotation(
        Placement(visible = true, transformation(origin = {108, 42}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.RealExpression off_freq(y = 0)  annotation(
        Placement(visible = true, transformation(origin = {74, 24}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.StateGraph.InitialStep initialStep(nIn=1, nOut= 1) annotation (Placement(visible = true, transformation(extent = {{-32, 132}, {-12, 152}}, rotation = 0)));
      Modelica.StateGraph.StepWithSignal step_on_off(nIn = 1, nOut = 1) annotation(
        Placement(visible = true, transformation(extent = {{24, 132}, {44, 152}}, rotation = 0)));
      Modelica.StateGraph.TransitionWithSignal transition2
        annotation (Placement(visible = true, transformation(extent = {{54, 132}, {74, 152}}, rotation = 0)));
        inner Modelica.StateGraph.StateGraphRoot stateGraphRoot
          annotation (Placement(visible = true, transformation(extent = {{-70, 136}, {-50, 156}}, rotation = 0)));
    Modelica.StateGraph.TransitionWithSignal transition1 annotation(
        Placement(visible = true, transformation(extent = {{-6, 132}, {14, 152}}, rotation = 0)));
    Modelica.Blocks.Interfaces.RealOutput y annotation(
        Placement(visible = true, transformation(origin = {148, 42}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {104, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Logical.GreaterThreshold greaterThreshold(threshold = dT)  annotation(
        Placement(visible = true, transformation(origin = {-22, 74}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
    Modelica.Blocks.Logical.Not not1 annotation(
        Placement(visible = true, transformation(origin = {20, 86}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Logical.And and1 annotation(
        Placement(visible = true, transformation(origin = {64, 112}, extent = {{-10, 10}, {10, -10}}, rotation = 90)));
    Modelica.Blocks.Logical.LessEqualThreshold lessEqualThreshold(threshold = lower)  annotation(
        Placement(visible = true, transformation(origin = {54, 72}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
    Modelica.Blocks.Nonlinear.Limiter limiterMin(uMax = 80, uMin = 20) annotation(
        Placement(visible = true, transformation(origin = {74, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(T_BT, feedback1.u2) annotation(
        Line(points = {{-104, -28}, {-40, -28}, {-40, 34}}, color = {0, 0, 127}));
      connect(T_set, feedback1.u1) annotation(
        Line(points = {{-102, 42}, {-48, 42}}, color = {0, 0, 127}));
      connect(feedback1.y, pi.u) annotation(
        Line(points = {{-30, 42}, {-14, 42}}, color = {0, 0, 127}));
    connect(pi.y, limiterMax.u) annotation(
        Line(points = {{10, 42}, {22, 42}}, color = {0, 0, 127}));
      connect(initialStep.outPort[1], transition1.inPort) annotation(
        Line(points = {{-11.5, 142}, {0.5, 142}}));
      connect(transition1.outPort, step_on_off.inPort[1]) annotation(
        Line(points = {{5.5, 142}, {23.5, 142}}));
      connect(step_on_off.outPort[1], transition2.inPort) annotation(
        Line(points = {{44.5, 142}, {60.5, 142}}));
      connect(transition2.outPort, initialStep.inPort[1]) annotation(
        Line(points = {{65.5, 142}, {83.5, 142}, {83.5, 160}, {-44.5, 160}, {-44.5, 142}, {-32.5, 142}}));
      connect(step_on_off.active, switch1.u2) annotation(
        Line(points = {{34, 131}, {34, 90}, {88, 90}, {88, 42}, {96, 42}}, color = {255, 0, 255}));
      connect(switch1.y, y) annotation(
        Line(points = {{119, 42}, {148, 42}}, color = {0, 0, 127}));
      connect(feedback1.y, greaterThreshold.u) annotation(
        Line(points = {{-30, 42}, {-22, 42}, {-22, 62}}, color = {0, 0, 127}));
      connect(greaterThreshold.y, transition1.condition) annotation(
        Line(points = {{-22, 86}, {4, 86}, {4, 130}}, color = {255, 0, 255}));
      connect(greaterThreshold.y, not1.u) annotation(
        Line(points = {{-22, 86}, {8, 86}}, color = {255, 0, 255}));
      connect(and1.y, transition2.condition) annotation(
        Line(points = {{64, 124}, {64, 130}}, color = {255, 0, 255}));
      connect(and1.u2, not1.y) annotation(
        Line(points = {{56, 100}, {56, 86}, {32, 86}}, color = {255, 0, 255}));
    connect(limiterMax.y, lessEqualThreshold.u) annotation(
        Line(points = {{46, 42}, {54, 42}, {54, 60}}, color = {0, 0, 127}));
    connect(lessEqualThreshold.y, and1.u1) annotation(
        Line(points = {{54, 84}, {64, 84}, {64, 100}}, color = {255, 0, 255}));
    connect(limiterMax.y, limiterMin.u) annotation(
        Line(points = {{46, 42}, {54, 42}, {54, 50}, {62, 50}}, color = {0, 0, 127}));
    connect(limiterMin.y, switch1.u1) annotation(
        Line(points = {{85, 50}, {96, 50}}, color = {0, 0, 127}));
    connect(off_freq.y, switch1.u3) annotation(
        Line(points = {{85, 24}, {90.5, 24}, {90.5, 34}, {96, 34}}, color = {0, 0, 127}));
      annotation(
        Icon(graphics = {Text(origin = {-48, -60}, extent = {{-32, 12}, {32, -12}}, textString = "T_BT"), Text(origin = {-42, 60}, extent = {{-34, 12}, {34, -12}}, textString = "T_set"), Rectangle(extent = {{-100, 100}, {100, -100}}), Text(origin = {72, -1}, extent = {{-18, 17}, {18, -17}}, textString = "f"), Text(origin = {-3, 120}, lineColor = {0, 0, 255}, extent = {{-89, 20}, {89, -20}}, textString = "%name")}),
        Diagram(coordinateSystem(extent = {{-120, 160}, {160, -60}})));
    end HexControl;
  end MyComponents;

  package Test
    model PlugFlowPipe
      extends Modelica.Icons.Example;
      replaceable package Medium = Buildings.Media.Water "Medium in the pipe" annotation(
        choicesAllMatching = true);
      final parameter Modelica.Units.SI.MassFlowRate m_flow_nominal = 6 "Mass flow rate";
      Modelica.Blocks.Sources.Ramp Tin(height = 20, duration = 0, offset = 273.15 + 50, startTime = 100) "Ramp pressure signal" annotation(
        Placement(transformation(extent = {{-100, -10}, {-80, 10}})));
      Buildings.Fluid.Sources.Boundary_pT sin(redeclare package Medium = Medium, T = 273.15 + 10, nPorts = 1, p(displayUnit = "Pa") = 101325) "Pressure boundary condition" annotation(
        Placement(transformation(extent = {{100, -10}, {80, 10}})));
      Buildings.Fluid.FixedResistances.PlugFlowPipe pip(redeclare package Medium = Medium, dh = 0.025, length = 0.5, m_flow_nominal = m_flow_nominal, cPip = 386, thickness = 0.005, initDelay = true, m_flow_start = m_flow_nominal, rhoPip = 8960, T_start_in = 323.15, T_start_out = 323.15, dIns = 0.005, kIns = 398) "Pipe" annotation(
        Placement(transformation(extent = {{0, 10}, {20, 30}})));
      Buildings.HeatTransfer.Sources.FixedTemperature bou[2](each T = 283.15) "Boundary temperature" annotation(
        Placement(transformation(extent = {{-40, 60}, {-20, 80}})));
      Buildings.Fluid.Sources.MassFlowSource_T sou(redeclare package Medium = Medium, use_T_in = true, m_flow = 0.17, nPorts = 1) "Flow source" annotation(
        Placement(transformation(extent = {{-60, 10}, {-40, 30}})));
      Buildings.Fluid.Sensors.TemperatureTwoPort senTemOut(redeclare package Medium = Medium, m_flow_nominal = m_flow_nominal, tau = 0, T_start = 323.15) "Temperature sensor" annotation(
        Placement(transformation(extent = {{40, 10}, {60, 30}})));
      Buildings.Fluid.Sensors.TemperatureTwoPort senTemIn(redeclare package Medium = Medium, m_flow_nominal = m_flow_nominal, tau = 0, T_start = 323.15) "Temperature sensor" annotation(
        Placement(transformation(extent = {{-30, 10}, {-10, 30}})));
      inner Modelica.Fluid.System system annotation(
        Placement(transformation(origin = {-90, 68}, extent = {{-10, -10}, {10, 10}})));
    equation
      connect(Tin.y, sou.T_in) annotation(
        Line(points = {{-79, 0}, {-68, 0}, {-68, 24}, {-62, 24}}, color = {0, 0, 127}));
      connect(pip.port_b, senTemOut.port_a) annotation(
        Line(points = {{20, 20}, {40, 20}}, color = {0, 127, 255}));
      connect(senTemOut.port_b, sin.ports[1]) annotation(
        Line(points = {{60, 20}, {76, 20}, {76, -1}, {80, -1}}, color = {0, 127, 255}));
      connect(senTemIn.port_b, pip.port_a) annotation(
        Line(points = {{-10, 20}, {0, 20}}, color = {0, 127, 255}));
      connect(bou[1].port, pip.heatPort) annotation(
        Line(points = {{-20, 70}, {-4, 70}, {-4, 40}, {10, 40}, {10, 30}}, color = {191, 0, 0}));
      connect(sou.ports[1], senTemIn.port_a) annotation(
        Line(points = {{-40, 20}, {-30, 20}}, color = {0, 127, 255}));
      annotation(
        experiment(StartTime = 0, StopTime = 1000, Tolerance = 1e-6, Interval = 2));
    end PlugFlowPipe;

    model PlugFlowPipeChangeP
      extends Modelica.Icons.Example;
      replaceable package Medium = Buildings.Media.Water "Medium in the pipe" annotation(
        choicesAllMatching = true);
      final parameter Modelica.Units.SI.MassFlowRate m_flow_nominal = 6 "Mass flow rate";
      Modelica.Blocks.Sources.Ramp Tin(height = 20, duration = 0, offset = 273.15 + 50, startTime = 100) "Ramp pressure signal" annotation(
        Placement(transformation(extent = {{-100, -10}, {-80, 10}})));
      Buildings.Fluid.Sources.Boundary_pT sin(redeclare package Medium = Medium, T = 273.15 + 10, nPorts = 2, p(displayUnit = "Pa") = 101325) "Pressure boundary condition" annotation(
        Placement(transformation(extent = {{100, -10}, {80, 10}})));
      Buildings.Fluid.FixedResistances.PlugFlowPipe pip(redeclare package Medium = Medium, dh = 0.1, length = 100, dIns = 0.05, kIns = 0.028, m_flow_nominal = m_flow_nominal, cPip = 500, thickness = 0.0032, initDelay = true, m_flow_start = m_flow_nominal, rhoPip = 8000, T_start_in = 323.15, T_start_out = 323.15, roughness = 2.5e-3) "Pipe" annotation(
        Placement(transformation(extent = {{0, 10}, {20, 30}})));
      Buildings.HeatTransfer.Sources.FixedTemperature bou[2](each T = 283.15) "Boundary temperature" annotation(
        Placement(transformation(extent = {{-40, 60}, {-20, 80}})));
      Buildings.Fluid.Sources.MassFlowSource_T sou(redeclare package Medium = Medium, use_T_in = true, m_flow = m_flow_nominal, nPorts = 1) "Flow source" annotation(
        Placement(transformation(extent = {{-60, 10}, {-40, 30}})));
      Buildings.Fluid.Sensors.TemperatureTwoPort senTemOut(redeclare package Medium = Medium, m_flow_nominal = m_flow_nominal, tau = 0, T_start = 323.15) "Temperature sensor" annotation(
        Placement(transformation(extent = {{40, 10}, {60, 30}})));
      Buildings.Fluid.Sensors.TemperatureTwoPort senTemIn(redeclare package Medium = Medium, m_flow_nominal = m_flow_nominal, tau = 0, T_start = 323.15) "Temperature sensor" annotation(
        Placement(transformation(extent = {{-30, 10}, {-10, 30}})));
    equation
      connect(Tin.y, sou.T_in) annotation(
        Line(points = {{-79, 0}, {-68, 0}, {-68, 24}, {-62, 24}}, color = {0, 0, 127}));
      connect(pip.port_b, senTemOut.port_a) annotation(
        Line(points = {{20, 20}, {40, 20}}, color = {0, 127, 255}));
      connect(senTemOut.port_b, sin.ports[1]) annotation(
        Line(points = {{60, 20}, {76, 20}, {76, -1}, {80, -1}}, color = {0, 127, 255}));
      connect(senTemIn.port_b, pip.port_a) annotation(
        Line(points = {{-10, 20}, {0, 20}}, color = {0, 127, 255}));
      connect(bou[1].port, pip.heatPort) annotation(
        Line(points = {{-20, 70}, {-4, 70}, {-4, 40}, {10, 40}, {10, 30}}, color = {191, 0, 0}));
      connect(sou.ports[1], senTemIn.port_a) annotation(
        Line(points = {{-40, 20}, {-30, 20}}, color = {0, 127, 255}));
    end PlugFlowPipeChangeP;

    model PlugFlowPipeChangePipe
      extends Modelica.Icons.Example;
      replaceable package Medium = Buildings.Media.Water "Medium in the pipe" annotation(
        choicesAllMatching = true);
      final parameter Modelica.Units.SI.MassFlowRate m_flow_nominal = 6 "Mass flow rate";
      Modelica.Blocks.Sources.Ramp Tin(height = 20, duration = 0, offset = 273.15 + 50, startTime = 100) "Ramp pressure signal" annotation(
        Placement(transformation(extent = {{-100, -10}, {-80, 10}})));
      Buildings.Fluid.Sources.Boundary_pT sin(redeclare package Medium = Medium, T = 273.15 + 10, nPorts = 1, p(displayUnit = "Pa") = 101325) "Pressure boundary condition" annotation(
        Placement(transformation(extent = {{100, -10}, {80, 10}})));
      Buildings.Fluid.FixedResistances.PlugFlowPipe pip(redeclare package Medium = Medium, dh = 0.025, length = 10, m_flow_nominal = m_flow_nominal, cPip = 386, thickness = 0.01, initDelay = true, m_flow_start = m_flow_nominal, rhoPip = 8960, T_start_in = 323.15, T_start_out = 323.15, dIns = 0.01, kIns = 100) "Pipe" annotation(
        Placement(transformation(extent = {{0, 10}, {20, 30}})));
      Buildings.HeatTransfer.Sources.FixedTemperature bou[2](each T = 283.15) "Boundary temperature" annotation(
        Placement(transformation(extent = {{-40, 60}, {-20, 80}})));
      Buildings.Fluid.Sources.MassFlowSource_T sou(redeclare package Medium = Medium, use_T_in = true, m_flow = 0.17, nPorts = 1) "Flow source" annotation(
        Placement(transformation(extent = {{-60, 10}, {-40, 30}})));
      Buildings.Fluid.Sensors.TemperatureTwoPort senTemOut(redeclare package Medium = Medium, m_flow_nominal = m_flow_nominal, tau = 0, T_start = 323.15) "Temperature sensor" annotation(
        Placement(transformation(extent = {{40, 10}, {60, 30}})));
      Buildings.Fluid.Sensors.TemperatureTwoPort senTemIn(redeclare package Medium = Medium, m_flow_nominal = m_flow_nominal, tau = 0, T_start = 323.15) "Temperature sensor" annotation(
        Placement(transformation(extent = {{-30, 10}, {-10, 30}})));
      inner Modelica.Fluid.System system annotation(
        Placement(transformation(origin = {-90, 68}, extent = {{-10, -10}, {10, 10}})));
    equation
      connect(Tin.y, sou.T_in) annotation(
        Line(points = {{-79, 0}, {-68, 0}, {-68, 24}, {-62, 24}}, color = {0, 0, 127}));
      connect(pip.port_b, senTemOut.port_a) annotation(
        Line(points = {{20, 20}, {40, 20}}, color = {0, 127, 255}));
      connect(senTemOut.port_b, sin.ports[1]) annotation(
        Line(points = {{60, 20}, {76, 20}, {76, -1}, {80, -1}}, color = {0, 127, 255}));
      connect(senTemIn.port_b, pip.port_a) annotation(
        Line(points = {{-10, 20}, {0, 20}}, color = {0, 127, 255}));
      connect(bou[1].port, pip.heatPort) annotation(
        Line(points = {{-20, 70}, {-4, 70}, {-4, 40}, {10, 40}, {10, 30}}, color = {191, 0, 0}));
      connect(sou.ports[1], senTemIn.port_a) annotation(
        Line(points = {{-40, 20}, {-30, 20}}, color = {0, 127, 255}));
      annotation(
        experiment(StartTime = 0, StopTime = 1000, Tolerance = 1e-06, Interval = 2));
    end PlugFlowPipeChangePipe;

    model Pump "Example model of movers using a parameter for setting the stage"
      extends Modelica.Icons.Example;
      package Medium = Buildings.Media.Water;
      parameter Modelica.Units.SI.MassFlowRate m_flow_nominal = 2 "Nominal mass flow rate";
      Modelica.Blocks.Sources.Ramp ramp(duration = 100, height = 0.9, offset = 0.1) "Ramp input for all movers" annotation(
        Placement(transformation(extent = {{-80, 60}, {-60, 80}})));
      Modelica.Blocks.Math.Gain gai_m_flow(k = m_flow_nominal) "Nominal mass flow rate" annotation(
        Placement(transformation(extent = {{-40, 10}, {-20, 30}})));
      Buildings.Fluid.Movers.FlowControlled_m_flow pump_m_flow(redeclare package Medium = Medium, m_flow_nominal = m_flow_nominal, use_inputFilter = false, massFlowRates = {0, 0.5, 1}*m_flow_nominal, inputType = Buildings.Fluid.Types.InputType.Continuous, energyDynamics = Modelica.Fluid.Types.Dynamics.SteadyState, nominalValuesDefineDefaultPressureCurve = true, allowFlowReversal = false) "Pump with m_flow input" annotation(
        Placement(transformation(extent = {{-10, -10}, {10, 10}})));
      Buildings.Fluid.Sources.Boundary_pT sou(redeclare package Medium = Medium, nPorts = 1) "Fluid source" annotation(
        Placement(transformation(extent = {{-80, -10}, {-60, 10}})));
      Buildings.Fluid.Sources.Boundary_pT sin(redeclare package Medium = Medium, nPorts = 1) "Fluid sink" annotation(
        Placement(transformation(extent = {{80, -10}, {60, 10}})));
      Buildings.Fluid.FixedResistances.PressureDrop res(redeclare package Medium = Medium, m_flow_nominal = m_flow_nominal, dp_nominal = dp_nominal) "Pressure drop component for avoiding singular system" annotation(
        Placement(transformation(origin = {2, 80}, extent = {{26, -90}, {46, -70}})));
      parameter Modelica.Units.SI.PressureDifference dp_nominal = 10000 "Nominal pressure raise";
    equation
      connect(gai_m_flow.u, ramp.y) annotation(
        Line(points = {{-42, 20}, {-50, 20}, {-50, 70}, {-59, 70}}, color = {0, 0, 127}));
      connect(sou.ports[1], pump_m_flow.port_a) annotation(
        Line(points = {{-60, -1.33333}, {-60, 0}, {-10, 0}}, color = {0, 127, 255}, smooth = Smooth.None));
      connect(res.port_b, sin.ports[1]) annotation(
        Line(points = {{48, 0}, {48, -80}, {60, -80}, {60, 1.33333}}, color = {0, 127, 255}));
      connect(pump_m_flow.port_b, res.port_a) annotation(
        Line(points = {{10, 0}, {28, 0}}, color = {0, 127, 255}));
      connect(gai_m_flow.y, pump_m_flow.m_flow_in) annotation(
        Line(points = {{-18, 20}, {0, 20}, {0, 12}}, color = {0, 0, 127}));
      annotation(
        Diagram);
    end Pump;

    model HEX_HP
      extends Modelica.Icons.Example;
      replaceable package Medium = Buildings.Media.Water "Medium in the pipe";
      inner Modelica.Fluid.System system annotation(
        Placement(transformation(origin = {-90, 70}, extent = {{-10, -10}, {10, 10}})));
      Modelica.Fluid.Sources.Boundary_pT boundary(redeclare package Medium = Medium, nPorts = 1, p = 101325, T = 273.15 + 10) annotation(
        Placement(transformation(origin = {82, 0}, extent = {{10, -10}, {-10, 10}})));
      Modelica.Fluid.Sources.MassFlowSource_T boundary1(redeclare package Medium = Medium, T = 273.15 + 10, nPorts = 1, use_m_flow_in = true) annotation(
        Placement(visible = true, transformation(origin = {-46, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Ramp ramp(duration = 10, height = 0.167, offset = 0, startTime = 0) annotation(
        Placement(visible = true, transformation(origin = {-86, 8}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Fluid.Sensors.Temperature temperature_a(redeclare package Medium = Medium) annotation(
        Placement(visible = true, transformation(origin = {-22, 18}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Fluid.Sensors.Temperature temperature_b(redeclare package Medium = Medium) annotation(
        Placement(transformation(origin = {38, 20}, extent = {{-10, -10}, {10, 10}})));
      MyComponents.HEX_HP hex_hp(K = 2, redeclare package Medium = Medium, T = 100, InitialTemp = 273.15 + 10) annotation(
        Placement(transformation(origin = {16, -22}, extent = {{-12, 14}, {12, 36}})));
      Modelica.Blocks.Sources.Ramp ramp1(duration = 0, height = 3000, offset = 0, startTime = 10) annotation(
        Placement(transformation(origin = {-32, 60}, extent = {{-10, -10}, {10, 10}})));
      Buildings.Fluid.FixedResistances.LosslessPipe pip1(redeclare package Medium = Medium, allowFlowReversal = false) annotation(
        Placement(transformation(origin = {52, 0}, extent = {{-10, -10}, {10, 10}})));
    equation
      connect(ramp.y, boundary1.m_flow_in) annotation(
        Line(points = {{-75, 8}, {-56, 8}}, color = {0, 0, 127}));
      connect(boundary1.ports[1], hex_hp.port_a) annotation(
        Line(points = {{-36, 0}, {6, 0}}, color = {0, 127, 255}));
      connect(hex_hp.port_a, temperature_a.port) annotation(
        Line(points = {{6, 0}, {-2, 0}, {-2, 8}, {-22, 8}}, color = {0, 127, 255}));
      connect(ramp1.y, hex_hp.freq) annotation(
        Line(points = {{-21, 60}, {-14, 60}, {-14, 9}, {6, 9}}, color = {0, 0, 127}));
      connect(temperature_b.port, hex_hp.port_b) annotation(
        Line(points = {{38, 10}, {38, 0}, {26, 0}}, color = {0, 127, 255}));
      connect(hex_hp.port_b, pip1.port_a) annotation(
        Line(points = {{26, 0}, {42, 0}}, color = {0, 127, 255}));
      connect(pip1.port_b, boundary.ports[1]) annotation(
        Line(points = {{62, 0}, {72, 0}}, color = {0, 127, 255}));
      annotation(
        experiment(StartTime = 0, StopTime = 150, Tolerance = 1e-6, Interval = 0.3),
        Diagram(coordinateSystem(extent = {{-100, 80}, {80, -20}})));
    end HEX_HP;

    model plant_test1
      extends Modelica.Icons.Example;
      replaceable package Medium = Buildings.Media.Water;
      inner Modelica.Fluid.System system annotation(
        Placement(transformation(origin = {-100, 134}, extent = {{-10, -10}, {10, 10}})));
      import SI = Modelica.Units.SI;
      // parameters
      parameter Modelica.Units.SI.Volume BT_vol = 0.001 "Buffer Tank Volume [m3]";
      parameter Modelica.Units.SI.Length BT_heiht = 0.1 "Buffer Tank height [m]";
      parameter Modelica.Units.SI.Temperature BT_initT = 273.15 + 30 "Buffer Tank Init Temperature [K]";
      parameter Modelica.Units.SI.Pressure BT_initP = 101325 "Buffer Tank Init Pressure [Pa]";
      parameter Modelica.Units.SI.Temperature Water_initT = 273.15 + 30;
      final parameter Modelica.Units.SI.MassFlowRate m_flow_nominal = 6 "Mass flow rate";
      // components
      Buildings.Fluid.Movers.FlowControlled_m_flow pump_m_flow(redeclare package Medium = Medium, m_flow_nominal = m_flow_nominal, use_inputFilter = false, massFlowRates = {0, 0.5, 1}*m_flow_nominal, inputType = Buildings.Fluid.Types.InputType.Continuous, energyDynamics = Modelica.Fluid.Types.Dynamics.SteadyState, nominalValuesDefineDefaultPressureCurve = true, allowFlowReversal = false) "Pump with m_flow input" annotation(
        Placement(transformation(origin = {-54, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
      MyComponents.HEX_HP hex_hp(InitialTemp = Water_initT, redeclare package Medium = Medium) "Heat Exchanger compresser frequency input" annotation(
        Placement(transformation(origin = {-74, 16}, extent = {{-12, -8}, {12, 14}}, rotation = 90)));
      Buildings.Fluid.FixedResistances.PlugFlowPipe pipe1(redeclare package Medium = Medium, length = 2, T_start_in = Water_initT, T_start_out = Water_initT, allowFlowReversal = false, cPip = 386, rhoPip = 8960, thickness = 0.005, dh = 0.025, dIns = 0.005, kIns = 100, m_flow_nominal = m_flow_nominal) "Outside PlugFlowPipe" annotation(
        Placement(transformation(origin = {-74, 46}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
      Buildings.Fluid.FixedResistances.PlugFlowPipe pipe4(redeclare package Medium = Medium, T_start_in = Water_initT, T_start_out = Water_initT, allowFlowReversal = false, cPip = 386, length = 2, rhoPip = 8960, thickness = 0.005, dh = 0.025, dIns = 0.005, kIns = 100, m_flow_nominal = m_flow_nominal) annotation(
        Placement(transformation(origin = {-74, -14}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
      Buildings.HeatTransfer.Sources.FixedTemperature outdoor(each T = 283.15) "Boundary temperature" annotation(
        Placement(transformation(origin = {-121, 113}, extent = {{-7, -7}, {7, 7}})));
      Modelica.Blocks.Sources.Ramp ramp(duration = 0, height = 3000, offset = 0, startTime = 10) annotation(
        Placement(transformation(origin = {-109, -5}, extent = {{-7, -7}, {7, 7}})));
      Modelica.Blocks.Sources.Ramp ramp1(duration = 30, height = 0.167, offset = 0, startTime = 0) annotation(
        Placement(transformation(origin = {-83, -61}, extent = {{-7, -7}, {7, 7}})));
      Buildings.Fluid.Sources.Boundary_ph bouOut(use_h_in = true, nPorts = 1, redeclare package Medium = Medium, use_X_in = false, use_Xi_in = false, use_C_in = false, use_p_in = false, p = 101325) annotation(
        Placement(transformation(origin = {16, 2}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      Buildings.Fluid.Sources.Boundary_pT bouIn(nPorts = 1, p = 101325, T = 303.15, redeclare package Medium = Medium, use_Xi_in = false, use_C_in = false, use_p_in = false, use_T_in = false) annotation(
        Placement(transformation(origin = {-8, 24}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
      Modelica.Fluid.Sensors.SpecificEnthalpyTwoPort specificEnthalpy(redeclare package Medium = Medium) annotation(
        Placement(transformation(origin = {-40, 66}, extent = {{-10, -10}, {10, 10}})));
    equation
      connect(pipe1.port_a, hex_hp.port_b) annotation(
        Line(points = {{-74, 36}, {-74, 26}}, color = {0, 127, 255}));
      connect(pipe4.port_b, hex_hp.port_a) annotation(
        Line(points = {{-74, -4}, {-74, 6}}, color = {0, 127, 255}));
      connect(pump_m_flow.port_b, pipe4.port_a) annotation(
        Line(points = {{-64, -40}, {-74, -40}, {-74, -24}}, color = {0, 127, 255}));
      connect(outdoor.port, pipe1.heatPort) annotation(
        Line(points = {{-114, 113}, {-90, 113}, {-90, 46}, {-84, 46}}, color = {191, 0, 0}));
      connect(outdoor.port, pipe4.heatPort) annotation(
        Line(points = {{-114, 113}, {-92, 113}, {-92, -14}, {-84, -14}}, color = {191, 0, 0}));
      connect(ramp.y, hex_hp.freq) annotation(
        Line(points = {{-101, -5}, {-82, -5}, {-82, 6}}, color = {0, 0, 127}));
      connect(ramp1.y, pump_m_flow.m_flow_in) annotation(
        Line(points = {{-76, -60}, {-54, -60}, {-54, -52}}, color = {0, 0, 127}));
      connect(bouOut.ports[1], pump_m_flow.port_a) annotation(
        Line(points = {{16, -8}, {16, -40}, {-44, -40}}, color = {0, 127, 255}));
      connect(pipe1.port_b, specificEnthalpy.port_a) annotation(
        Line(points = {{-74, 56}, {-74, 66}, {-50, 66}}, color = {0, 127, 255}));
      connect(specificEnthalpy.port_b, bouIn.ports[1]) annotation(
        Line(points = {{-30, 66}, {-8, 66}, {-8, 34}}, color = {0, 127, 255}));
      connect(specificEnthalpy.h_out, bouOut.h_in) annotation(
        Line(points = {{-40, 78}, {20, 78}, {20, 14}}, color = {0, 0, 127}));
      annotation(
        Diagram(coordinateSystem(extent = {{-140, 160}, {40, -80}})));
    end plant_test1;

    model plant_test2
      extends Modelica.Icons.Example;
      replaceable package Medium = Buildings.Media.Water;
      inner Modelica.Fluid.System system annotation(
        Placement(transformation(origin = {-124, 144}, extent = {{-10, -10}, {10, 10}})));
      import SI = Modelica.Units.SI;
      // Parameters
      // Initialization
      parameter Modelica.Units.SI.Temperature medium_initT = 273.15 + 30 "Water Initial Temperature [K]";
      parameter Modelica.Units.SI.Temperature BT_initT = 273.15 + 30 "Buffer Tank Init Temperature [K]";
      parameter Modelica.Units.SI.Pressure BT_initP = 101325 "Buffer Tank Init Pressure [Pa]";
      // Buffer Tank
      parameter Modelica.Units.SI.Volume BT_vol = 0.001 "Buffer Tank Volume [m3]";
      parameter Modelica.Units.SI.Length BT_height = 0.1 "Buffer Tank height [m]";
      // Pipe
      parameter Modelica.Units.SI.SpecificHeatCapacity pip_c = 386 "Specific heat of pipe wall material [J/kg.K]";
      parameter Modelica.Units.SI.Density pip_rho = 8960 "Density of pipe wall material [kg/m3]";
      parameter Modelica.Units.SI.Length pip_dh = 0.025 "Inner diameter of pipe [m]";
      parameter Modelica.Units.SI.Length pip_thickness = 0.005 "Thickness of pipe wall [m]";
      parameter Modelica.Units.SI.Length pip_dIns = 0.005 "Thickness of pipe insulation [m]";
      parameter Modelica.Units.SI.ThermalConductivity pip_kIns = 1 "Heat conductivity of pipe insulation [W/m.K]";
      //
      final parameter Modelica.Units.SI.MassFlowRate m_flow_nominal = 6 "Mass flow rate [kg/s]";
      final parameter Modelica.Units.SI.PressureDifference dp_nominal = 4500 "Design pressure drop";
      //===========
      // components
      Buildings.Fluid.Movers.FlowControlled_m_flow pump_m_flow(redeclare package Medium = Medium, m_flow_nominal = m_flow_nominal, use_inputFilter = false, massFlowRates = {0, 0.5, 1}*m_flow_nominal, inputType = Buildings.Fluid.Types.InputType.Continuous, energyDynamics = Modelica.Fluid.Types.Dynamics.SteadyState, nominalValuesDefineDefaultPressureCurve = true, allowFlowReversal = false, T_start = medium_initT) "Pump with m_flow input" annotation(
        Placement(transformation(origin = {-54, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
      MyComponents.HEX_HP hex_hp(InitialTemp = medium_initT, redeclare package Medium = Medium) "Heat Exchanger compresser frequency input" annotation(
        Placement(transformation(origin = {-74, 16}, extent = {{-12, -8}, {12, 14}}, rotation = 90)));
      Buildings.Fluid.FixedResistances.PlugFlowPipe pipe1(redeclare package Medium = Medium, length = 2, T_start_in = medium_initT, T_start_out = medium_initT, allowFlowReversal = false, cPip = pip_c, rhoPip = pip_rho, thickness = pip_thickness, dh = pip_dh, dIns = pip_dIns, kIns = pip_kIns, m_flow_nominal = m_flow_nominal) "Outside PlugFlowPipe" annotation(
        Placement(transformation(origin = {-74, 46}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
      Buildings.Fluid.Storage.Stratified BufferTank(redeclare package Medium = Medium, T_start = BT_initT, VTan = BT_vol, dIns = 0.01, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, hTan = BT_height, nSeg = 5, p_start = BT_initP) "comment" annotation(
        Placement(transformation(origin = {21, 77}, extent = {{15, -15}, {-15, 15}}, rotation = -180)));
      Buildings.Fluid.FixedResistances.PlugFlowPipe pipe3(redeclare package Medium = Medium, T_start_in = medium_initT, T_start_out = medium_initT, allowFlowReversal = false, cPip = pip_c, length = 5, rhoPip = pip_rho, thickness = pip_thickness, dh = pip_dh, dIns = pip_dIns, kIns = pip_kIns, m_flow_nominal = m_flow_nominal) "comment" annotation(
        Placement(transformation(origin = {26, 40}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      Buildings.Fluid.FixedResistances.PlugFlowPipe pipe2(redeclare package Medium = Medium, T_start_in = medium_initT, T_start_out = medium_initT, allowFlowReversal = false, cPip = pip_c, length = 2, rhoPip = pip_rho, thickness = pip_thickness, dh = pip_dh, dIns = pip_dIns, kIns = pip_kIns, m_flow_nominal = m_flow_nominal) "Indoor PlugFlowPipe" annotation(
        Placement(transformation(origin = {-24, 96}, extent = {{-10, -10}, {10, 10}})));
      Buildings.Fluid.FixedResistances.PlugFlowPipe pipe4(redeclare package Medium = Medium, T_start_in = medium_initT, T_start_out = medium_initT, allowFlowReversal = false, cPip = pip_c, length = 2, rhoPip = pip_rho, thickness = pip_thickness, dh = pip_dh, dIns = pip_dIns, kIns = pip_kIns, m_flow_nominal = m_flow_nominal) "" annotation(
        Placement(transformation(origin = {-74, -14}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
      Buildings.Fluid.Actuators.Valves.TwoWayLinear val(redeclare package Medium = Medium, m_flow_nominal = m_flow_nominal, CvData = Buildings.Fluid.Types.CvTypes.OpPoint, dpValve_nominal(displayUnit = "kPa") = dp_nominal, use_inputFilter = false) "Valve model, linear opening characteristics" annotation(
        Placement(transformation(origin = {-48, 66}, extent = {{-10, -10}, {10, 10}})));
      Buildings.Fluid.FixedResistances.PlugFlowPipe pipe5(redeclare package Medium = Medium, T_start_in = medium_initT, T_start_out = medium_initT, allowFlowReversal = false, cPip = pip_c, length = 0.2, rhoPip = pip_rho, thickness = pip_thickness, dh = pip_dh, dIns = pip_dIns, kIns = pip_kIns, m_flow_nominal = m_flow_nominal) "" annotation(
        Placement(transformation(origin = {-24, 42}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      Modelica.Blocks.Sources.Constant const(k = 0) annotation(
        Placement(transformation(origin = {-62, 86}, extent = {{-6, -6}, {6, 6}})));
      Buildings.HeatTransfer.Sources.FixedTemperature outdoor(each T = 283.15) "Boundary temperature" annotation(
        Placement(transformation(origin = {-121, 113}, extent = {{-7, -7}, {7, 7}})));
      Modelica.Blocks.Sources.Ramp ramp(duration = 0, height = 3000, offset = 0, startTime = 10) annotation(
        Placement(transformation(origin = {-109, -5}, extent = {{-7, -7}, {7, 7}})));
      Modelica.Blocks.Sources.Ramp ramp1(duration = 1, height = 0.167, offset = 0, startTime = 0) annotation(
        Placement(transformation(origin = {-85, -61}, extent = {{-7, -7}, {7, 7}})));
      Modelica.Fluid.Sensors.SpecificEnthalpyTwoPort specificEnthalpy(redeclare package Medium = Medium) annotation(
        Placement(transformation(origin = {26, 12}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      Buildings.Fluid.Storage.ExpansionVessel exp(redeclare package Medium = Medium, p_start = 101325) annotation(
        Placement(transformation(origin = {48, -20}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    equation
      connect(pipe1.port_a, hex_hp.port_b) annotation(
        Line(points = {{-74, 36}, {-74, 26}}, color = {0, 127, 255}));
      connect(pipe1.port_b, pipe2.port_a) annotation(
        Line(points = {{-74, 56}, {-74, 96}, {-34, 96}}, color = {0, 127, 255}));
      connect(pipe4.port_b, hex_hp.port_a) annotation(
        Line(points = {{-74, -4}, {-74, 6}}, color = {0, 127, 255}));
      connect(pump_m_flow.port_b, pipe4.port_a) annotation(
        Line(points = {{-64, -40}, {-74, -40}, {-74, -24}}, color = {0, 127, 255}));
      connect(pipe1.port_b, val.port_a) annotation(
        Line(points = {{-74, 56}, {-74, 66}, {-58, 66}}, color = {0, 127, 255}));
      connect(val.port_b, pipe5.port_a) annotation(
        Line(points = {{-38, 66}, {-24, 66}, {-24, 52}}, color = {0, 127, 255}));
      connect(pipe5.port_b, pump_m_flow.port_a) annotation(
        Line(points = {{-24, 32}, {-24, -40}, {-44, -40}}, color = {0, 127, 255}));
      connect(const.y, val.y) annotation(
        Line(points = {{-56, 86}, {-48, 86}, {-48, 78}}, color = {0, 0, 127}));
      connect(outdoor.port, pipe2.heatPort) annotation(
        Line(points = {{-114, 113}, {-24, 113}, {-24, 106}}, color = {191, 0, 0}));
      connect(outdoor.port, pipe1.heatPort) annotation(
        Line(points = {{-114, 113}, {-90, 113}, {-90, 46}, {-84, 46}}, color = {191, 0, 0}));
      connect(outdoor.port, pipe4.heatPort) annotation(
        Line(points = {{-114, 113}, {-92, 113}, {-92, -14}, {-84, -14}}, color = {191, 0, 0}));
      connect(outdoor.port, pipe5.heatPort) annotation(
        Line(points = {{-114, 113}, {-6, 113}, {-6, 42}, {-14, 42}}, color = {191, 0, 0}));
      connect(outdoor.port, pipe3.heatPort) annotation(
        Line(points = {{-114, 113}, {46, 113}, {46, 40}, {36, 40}}, color = {191, 0, 0}));
      connect(ramp.y, hex_hp.freq) annotation(
        Line(points = {{-101, -5}, {-82, -5}, {-82, 6}}, color = {0, 0, 127}));
      connect(ramp1.y, pump_m_flow.m_flow_in) annotation(
        Line(points = {{-77, -61}, {-77, -52}, {-54, -52}}, color = {0, 0, 127}));
      connect(pipe3.port_b, specificEnthalpy.port_a) annotation(
        Line(points = {{26, 30}, {26, 22}}, color = {0, 127, 255}));
      connect(pipe2.port_b, BufferTank.port_a) annotation(
        Line(points = {{-14, 96}, {6, 96}, {6, 77}}, color = {0, 127, 255}));
      connect(BufferTank.port_b, pipe3.port_a) annotation(
        Line(points = {{36, 77}, {42, 77}, {42, 56}, {26, 56}, {26, 50}}, color = {0, 127, 255}));
      connect(exp.port_a, specificEnthalpy.port_b) annotation(
        Line(points = {{38, -20}, {26, -20}, {26, 2}}, color = {0, 127, 255}));
      connect(specificEnthalpy.port_b, pump_m_flow.port_a) annotation(
        Line(points = {{26, 2}, {26, -40}, {-44, -40}}, color = {0, 127, 255}));
      annotation(
        Diagram(coordinateSystem(extent = {{-140, 160}, {60, -80}})),
        experiment(StartTime = 0, StopTime = 10000, Tolerance = 1e-6, Interval = 20));
    end plant_test2;

    model Room01 "Example model of movers using a parameter for setting the stage"
      extends Modelica.Icons.Example;
      //Parameter
      inner Modelica.Fluid.System system annotation(
        Placement(transformation(origin = {-56, 34}, extent = {{-10, -10}, {10, 10}})));
      package Medium = Buildings.Media.Water;
      parameter Modelica.Units.SI.PressureDifference dp_nominal = 10000 "Nominal pressure raise";
      //Initialize
      parameter Modelica.Units.SI.Temperature medium_initT = 273.15 + 30 "Water Initial Temperature [K]";
      parameter Modelica.Units.SI.Temperature rad_initT = 273.15 + 40;
      //
      parameter Modelica.Units.SI.SpecificHeatCapacity pip_c = 386 "Specific heat of pipe wall material [J/kg.K]";
      parameter Modelica.Units.SI.Density pip_rho = 8960 "Density of pipe wall material [kg/m3]";
      parameter Modelica.Units.SI.Length pip_dh = 0.025 "Inner diameter of pipe [m]";
      parameter Modelica.Units.SI.Length pip_thickness = 0.005 "Thickness of pipe wall [m]";
      parameter Modelica.Units.SI.Length pip_dIns = 0.005 "Thickness of pipe insulation [m]";
      parameter Modelica.Units.SI.ThermalConductivity pip_kIns = 1 "Heat conductivity of pipe insulation [W/m.K]";
      //Components
      parameter Modelica.Units.SI.MassFlowRate m_flow_nominal = 2 "Nominal mass flow rate";
      Buildings.Fluid.Sources.Boundary_pT sou(redeclare package Medium = Medium, nPorts = 1, p(displayUnit = "Pa") = 101325, T = 273.15 + 50) "Fluid source" annotation(
        Placement(transformation(origin = {2, 10}, extent = {{-80, -10}, {-60, 10}})));
      MyComponents.roomCycle roomCycle1 annotation(
        Placement(transformation(origin = {-6, 0}, extent = {{-20, -20}, {20, 20}})));
      Buildings.Fluid.Sources.Boundary_pT sin(redeclare package Medium = Medium, nPorts = 1, p(displayUnit = "Pa") = 101325, T = 273.15 + 30) "Fluid sink" annotation(
        Placement(transformation(origin = {2, -20}, extent = {{-80, -10}, {-60, 10}})));
    equation
      connect(sou.ports[1], roomCycle1.port_a) annotation(
        Line(points = {{-58, 10}, {-26, 10}}, color = {0, 127, 255}));
      connect(sin.ports[1], roomCycle1.port_b) annotation(
        Line(points = {{-58, -20}, {-26, -20}, {-26, -10}}, color = {0, 127, 255}));
      annotation(
        Diagram(coordinateSystem(extent = {{-80, 60}, {20, -40}})),
        experiment(StartTime = 0, StopTime = 100, Tolerance = 1e-6, Interval = 0.2));
    end Room01;

    model Room02 "Example model of movers using a parameter for setting the stage"
      extends Modelica.Icons.Example;
      //Parameter
      inner Modelica.Fluid.System system annotation(
        Placement(transformation(origin = {-56, 34}, extent = {{-10, -10}, {10, 10}})));
      package Medium = Buildings.Media.Water;
      parameter Modelica.Units.SI.PressureDifference dp_nominal = 10000 "Nominal pressure raise";
      //Initialize
      parameter Modelica.Units.SI.Temperature medium_initT = 273.15 + 30 "Water Initial Temperature [K]";
      parameter Modelica.Units.SI.Temperature rad_initT = 273.15 + 40;
      //
      parameter Modelica.Units.SI.SpecificHeatCapacity pip_c = 386 "Specific heat of pipe wall material [J/kg.K]";
      parameter Modelica.Units.SI.Density pip_rho = 8960 "Density of pipe wall material [kg/m3]";
      parameter Modelica.Units.SI.Length pip_dh = 0.025 "Inner diameter of pipe [m]";
      parameter Modelica.Units.SI.Length pip_thickness = 0.005 "Thickness of pipe wall [m]";
      parameter Modelica.Units.SI.Length pip_dIns = 0.005 "Thickness of pipe insulation [m]";
      parameter Modelica.Units.SI.ThermalConductivity pip_kIns = 1 "Heat conductivity of pipe insulation [W/m.K]";
      //Components
      parameter Modelica.Units.SI.MassFlowRate m_flow_nominal = 2 "Nominal mass flow rate";
      Buildings.Fluid.Sources.Boundary_pT sou(redeclare package Medium = Medium, nPorts = 1, p(displayUnit = "Pa") = 101325, T = 273.15 + 50) "Fluid source" annotation(
        Placement(transformation(origin = {2, 10}, extent = {{-80, -10}, {-60, 10}})));
      Buildings.Fluid.Sources.Boundary_pT sin(redeclare package Medium = Medium, nPorts = 1, p(displayUnit = "Pa") = 101325, T = 273.15 + 30) "Fluid sink" annotation(
        Placement(transformation(origin = {2, -20}, extent = {{-80, -10}, {-60, 10}})));
      MyComponents.roomCycleWithoutLossnay roomCycleWithoutLossnay1 annotation(
        Placement(transformation(origin = {-14, -6}, extent = {{-10, -10}, {10, 10}})));
    equation
      connect(sou.ports[1], roomCycleWithoutLossnay1.port_a) annotation(
        Line(points = {{-58, 10}, {-24, 10}, {-24, 0}}, color = {0, 127, 255}));
      connect(roomCycleWithoutLossnay1.port_b, sin.ports[1]) annotation(
        Line(points = {{-24, -10}, {-58, -10}, {-58, -20}}, color = {0, 127, 255}));
      annotation(
        Diagram(coordinateSystem(extent = {{-80, 60}, {20, -40}})),
        experiment(StartTime = 0, StopTime = 100, Tolerance = 1e-6, Interval = 0.2));
    end Room02;

    model Lossnay "Example model of movers using a parameter for setting the stage"
      extends Modelica.Icons.Example;
      replaceable package MediumAir = Buildings.Media.Air;
      //Parameter
      inner Modelica.Fluid.System system annotation(
        Placement(transformation(origin = {-56, 34}, extent = {{-10, -10}, {10, 10}})));
      package Medium = Buildings.Media.Water;
      parameter Modelica.Units.SI.PressureDifference dp_nominal = 10000 "Nominal pressure raise";
      parameter Modelica.Units.SI.Temperature mediumRoomAir_initT = 273.15 + 30 "RoomAir Initial Temperature [K]";
      //Initialize
      parameter Modelica.Units.SI.Temperature medium_initT = 273.15 + 30 "Water Initial Temperature [K]";
      parameter Modelica.Units.SI.Temperature rad_initT = 273.15 + 40;
      //
      parameter Modelica.Units.SI.SpecificHeatCapacity pip_c = 386 "Specific heat of pipe wall material [J/kg.K]";
      parameter Modelica.Units.SI.Density pip_rho = 8960 "Density of pipe wall material [kg/m3]";
      parameter Modelica.Units.SI.Length pip_dh = 0.025 "Inner diameter of pipe [m]";
      parameter Modelica.Units.SI.Length pip_thickness = 0.005 "Thickness of pipe wall [m]";
      parameter Modelica.Units.SI.Length pip_dIns = 0.005 "Thickness of pipe insulation [m]";
      parameter Modelica.Units.SI.ThermalConductivity pip_kIns = 1 "Heat conductivity of pipe insulation [W/m.K]";
      //Components
      parameter Modelica.Units.SI.MassFlowRate m_flow_nominal = 2 "Nominal mass flow rate";
      MyComponents.Lossnay lossnay annotation(
        Placement(transformation(origin = {-4, 2}, extent = {{-10, -10}, {10, 10}})));
      Buildings.Fluid.MixingVolumes.MixingVolume roomAir(redeclare package Medium = MediumAir, m_flow_nominal = 0.01, V = 10, p_start(displayUnit = "Pa") = 102315, T_start = 273.15 + 20, nPorts = 2, massDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial) annotation(
        Placement(transformation(origin = {-52, 6}, extent = {{-10, -10}, {10, 10}})));
    equation
      connect(roomAir.ports[1], lossnay.port_a) annotation(
        Line(points = {{-52, -4}, {-52, -6}, {-20, -6}, {-20, 8}, {-14, 8}}, color = {0, 127, 255}));
      connect(roomAir.ports[2], lossnay.port_b) annotation(
        Line(points = {{-52, -4}, {-52, -12}, {-16, -12}, {-16, -4}, {-14, -4}}, color = {0, 127, 255}));
      annotation(
        Diagram(coordinateSystem(extent = {{-80, 60}, {20, -40}})),
        experiment(StartTime = 0, StopTime = 100, Tolerance = 1e-6, Interval = 0.2));
    end Lossnay;
    
    model FirstExample_Variant2
        "A variant of the first simple StateGraph example"
      extends Modelica.Icons.Example;
      Modelica.StateGraph.InitialStep initialStep(nIn=1, nOut=1) annotation (Placement(transformation(extent={{-70,0},{-50,20}})));
      Modelica.StateGraph.Transition transition1(condition = false,enableTimer= true, waitTime= 0.1)
        annotation (Placement(transformation(extent={{-42,0},{-22,20}})));
      Modelica.StateGraph.StepWithSignal step(nIn=1, nOut=1)
                annotation (Placement(transformation(extent={{-14,0},{6,20}})));
      Modelica.StateGraph.TransitionWithSignal transition2
        annotation (Placement(transformation(extent={{52,0},{72,20}})));
      Modelica.Blocks.Logical.Timer timer annotation (Placement(transformation(
                extent={{6,-40},{26,-20}})));
      Modelica.Blocks.Logical.GreaterEqualThreshold greaterEqual(threshold=1)
        annotation (Placement(transformation(extent={{36,-40},{56,-20}})));
        inner Modelica.StateGraph.StateGraphRoot stateGraphRoot
          annotation (Placement(transformation(extent={{-80,60},{-60,80}})));
    equation
    
      connect(initialStep.outPort[1], transition1.inPort)
        annotation (Line(points={{-49.5,10},{-36,10}}));
    
      connect(transition1.outPort, step.inPort[1])
        annotation (Line(points={{-30.5,10},{-15,10}}));
      connect(step.active, timer.u) annotation (Line(points={{-4,-1},{-4,-30},{4,
                -30}}, color={255,0,255}));
      connect(step.outPort[1], transition2.inPort)
        annotation (Line(points={{6.5,10},{58,10}}));
      connect(timer.y, greaterEqual.u)
        annotation (Line(points={{27,-30},{34,-30}}, color={0,0,255}));
      connect(greaterEqual.y, transition2.condition) annotation (Line(points={{57,
                -30},{62,-30},{62,-2}}, color={255,0,255}));
      connect(transition2.outPort, initialStep.inPort[1]) annotation (Line(points=
               {{63.5,10},{82,10},{82,32},{-80,32},{-80,10},{-71,10}}));
      annotation (experiment(StopTime=5.5));
    end FirstExample_Variant2;
    
    model FirstExample_Variant2_edit1
        "A variant of the first simple StateGraph example"
      extends Modelica.Icons.Example;
      Modelica.StateGraph.InitialStep initialStep(nIn=1, nOut= 1) annotation (Placement(transformation(extent={{-70,0},{-50,20}})));
      Modelica.StateGraph.StepWithSignal step(nIn = 1, nOut = 1) annotation(
        Placement(transformation(extent = {{-14, 0}, {6, 20}})));
      Modelica.StateGraph.TransitionWithSignal transition2
        annotation (Placement(visible = true, transformation(extent = {{16, 0}, {36, 20}}, rotation = 0)));
      Modelica.Blocks.Logical.Timer timer annotation (Placement(transformation(
                extent={{6,-40},{26,-20}})));
      Modelica.Blocks.Logical.GreaterEqualThreshold greaterEqual(threshold=1)
        annotation (Placement(transformation(extent={{36,-40},{56,-20}})));
        inner Modelica.StateGraph.StateGraphRoot stateGraphRoot
          annotation (Placement(transformation(extent={{-80,60},{-60,80}})));
  Modelica.StateGraph.TransitionWithSignal transition1 annotation(
        Placement(visible = true, transformation(extent = {{-44, 0}, {-24, 20}}, rotation = 0)));
  Modelica.Blocks.Sources.BooleanPulse booleanPulse(period = 2, startTime = 1, width = 10)  annotation(
        Placement(visible = true, transformation(origin = {-70, -26}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(step.active, timer.u) annotation(
        Line(points = {{-4, -1}, {-4, -30}, {4, -30}}, color = {255, 0, 255}));
  connect(step.outPort[1], transition2.inPort) annotation(
        Line(points = {{6.5, 10}, {22, 10}}));
      connect(timer.y, greaterEqual.u) annotation(
        Line(points = {{27, -30}, {34, -30}}, color = {0, 0, 255}));
  connect(greaterEqual.y, transition2.condition) annotation(
        Line(points = {{57, -30}, {26, -30}, {26, -2}}, color = {255, 0, 255}));
  connect(transition2.outPort, initialStep.inPort[1]) annotation(
        Line(points = {{27.5, 10}, {40, 10}, {40, 32}, {-80, 32}, {-80, 10}, {-71, 10}}));
  connect(initialStep.outPort[1], transition1.inPort) annotation(
        Line(points = {{-50, 10}, {-38, 10}}));
  connect(transition1.outPort, step.inPort[1]) annotation(
        Line(points = {{-32, 10}, {-14, 10}}));
  connect(booleanPulse.y, transition1.condition) annotation(
        Line(points = {{-58, -26}, {-34, -26}, {-34, -2}}, color = {255, 0, 255}));
      annotation (experiment(StopTime=5.5));
    end FirstExample_Variant2_edit1;
    
    model HexControl
      extends Modelica.Icons.Example;
    System.MyComponents.HexControl hexControl annotation(
          Placement(visible = true, transformation(origin = {-6, 14}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.RealExpression realExpression(y = 30)  annotation(
          Placement(visible = true, transformation(origin = {-48, 22}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Sine sine(amplitude = 35, f = 0.01, startTime = 1)  annotation(
          Placement(visible = true, transformation(origin = {-48, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(realExpression.y, hexControl.T_set) annotation(
        Line(points = {{-36, 22}, {-16, 22}, {-16, 20}}, color = {0, 0, 127}));
  connect(sine.y, hexControl.T_BT) annotation(
        Line(points = {{-36, 2}, {-22, 2}, {-22, 8}, {-16, 8}}, color = {0, 0, 127}));
      annotation (experiment(StopTime = 100, StartTime = 0, Tolerance = 1e-06, Interval = 0.2));
    end HexControl;
  end Test;

  model plantold1
    replaceable package Medium = Buildings.Media.Water;
    inner Modelica.Fluid.System system annotation(
      Placement(transformation(origin = {-124, 144}, extent = {{-10, -10}, {10, 10}})));
    import SI = Modelica.Units.SI;
    // Parameters
    // Initialization
    parameter Modelica.Units.SI.Temperature medium_initT = 273.15 + 30 "Water Initial Temperature [K]";
    parameter Modelica.Units.SI.Temperature BT_initT = 273.15 + 30 "Buffer Tank Init Temperature [K]";
    parameter Modelica.Units.SI.Pressure BT_initP = 101325 "Buffer Tank Init Pressure [Pa]";
    // Buffer Tank
    parameter Modelica.Units.SI.Volume BT_vol = 0.001 "Buffer Tank Volume [m3]";
    parameter Modelica.Units.SI.Length BT_height = 0.1 "Buffer Tank height [m]";
    // Pipe
    parameter Modelica.Units.SI.SpecificHeatCapacity pip_c = 386 "Specific heat of pipe wall material [J/kg.K]";
    parameter Modelica.Units.SI.Density pip_rho = 8960 "Density of pipe wall material [kg/m3]";
    parameter Modelica.Units.SI.Length pip_dh = 0.025 "Inner diameter of pipe [m]";
    parameter Modelica.Units.SI.Length pip_thickness = 0.005 "Thickness of pipe wall [m]";
    parameter Modelica.Units.SI.Length pip_dIns = 0.005 "Thickness of pipe insulation [m]";
    parameter Modelica.Units.SI.ThermalConductivity pip_kIns = 1 "Heat conductivity of pipe insulation [W/m.K]";
    //
    final parameter Modelica.Units.SI.MassFlowRate m_flow_nominal = 6 "Mass flow rate [kg/s]";
    final parameter Modelica.Units.SI.PressureDifference dp_nominal = 4500 "Design pressure drop";
    //===========
    // components
    Buildings.Fluid.Movers.FlowControlled_m_flow pump_m_flow(redeclare package Medium = Medium, m_flow_nominal = m_flow_nominal, use_inputFilter = false, massFlowRates = {0, 0.5, 1}*m_flow_nominal, inputType = Buildings.Fluid.Types.InputType.Continuous, energyDynamics = Modelica.Fluid.Types.Dynamics.SteadyState, nominalValuesDefineDefaultPressureCurve = true, allowFlowReversal = false, T_start = medium_initT) "Pump with m_flow input" annotation(
      Placement(transformation(origin = {-54, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
    MyComponents.HEX_HP hex_hp(InitialTemp = medium_initT, redeclare package Medium = Medium) "Heat Exchanger compresser frequency input" annotation(
      Placement(transformation(origin = {-74, 16}, extent = {{-12, -8}, {12, 14}}, rotation = 90)));
    Buildings.Fluid.FixedResistances.PlugFlowPipe pipe1(redeclare package Medium = Medium, length = 2, T_start_in = medium_initT, T_start_out = medium_initT, allowFlowReversal = false, cPip = pip_c, rhoPip = pip_rho, thickness = pip_thickness, dh = pip_dh, dIns = pip_dIns, kIns = pip_kIns, m_flow_nominal = m_flow_nominal) "Outside PlugFlowPipe" annotation(
      Placement(transformation(origin = {-74, 46}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
    Buildings.Fluid.Storage.Stratified BufferTank(redeclare package Medium = Medium, T_start = BT_initT, VTan = BT_vol, dIns = 0.01, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, hTan = BT_height, nSeg = 5, p_start = BT_initP) "comment" annotation(
      Placement(transformation(origin = {21, 77}, extent = {{15, -15}, {-15, 15}}, rotation = -180)));
    Buildings.Fluid.FixedResistances.PlugFlowPipe pipe3(redeclare package Medium = Medium, T_start_in = medium_initT, T_start_out = medium_initT, allowFlowReversal = false, cPip = pip_c, length = 5, rhoPip = pip_rho, thickness = pip_thickness, dh = pip_dh, dIns = pip_dIns, kIns = pip_kIns, m_flow_nominal = m_flow_nominal) "comment" annotation(
      Placement(transformation(origin = {26, 40}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Buildings.Fluid.FixedResistances.PlugFlowPipe pipe2(redeclare package Medium = Medium, T_start_in = medium_initT, T_start_out = medium_initT, allowFlowReversal = false, cPip = pip_c, length = 2, rhoPip = pip_rho, thickness = pip_thickness, dh = pip_dh, dIns = pip_dIns, kIns = pip_kIns, m_flow_nominal = m_flow_nominal) "Indoor PlugFlowPipe" annotation(
      Placement(transformation(origin = {-24, 96}, extent = {{-10, -10}, {10, 10}})));
    Buildings.Fluid.FixedResistances.PlugFlowPipe pipe4(redeclare package Medium = Medium, T_start_in = medium_initT, T_start_out = medium_initT, allowFlowReversal = false, cPip = pip_c, length = 2, rhoPip = pip_rho, thickness = pip_thickness, dh = pip_dh, dIns = pip_dIns, kIns = pip_kIns, m_flow_nominal = m_flow_nominal) "" annotation(
      Placement(transformation(origin = {-74, -14}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
    Buildings.Fluid.Actuators.Valves.TwoWayLinear val(redeclare package Medium = Medium, m_flow_nominal = m_flow_nominal, CvData = Buildings.Fluid.Types.CvTypes.OpPoint, dpValve_nominal(displayUnit = "kPa") = dp_nominal, use_inputFilter = false) "Valve model, linear opening characteristics" annotation(
      Placement(transformation(origin = {-48, 66}, extent = {{-10, -10}, {10, 10}})));
    Buildings.Fluid.FixedResistances.PlugFlowPipe pipe5(redeclare package Medium = Medium, T_start_in = medium_initT, T_start_out = medium_initT, allowFlowReversal = false, cPip = pip_c, length = 0.2, rhoPip = pip_rho, thickness = pip_thickness, dh = pip_dh, dIns = pip_dIns, kIns = pip_kIns, m_flow_nominal = m_flow_nominal) "" annotation(
      Placement(transformation(origin = {-24, 42}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Blocks.Sources.Constant const(k = 0) annotation(
      Placement(transformation(origin = {-62, 86}, extent = {{-6, -6}, {6, 6}})));
    Buildings.HeatTransfer.Sources.FixedTemperature outdoor(each T = 283.15) "Boundary temperature" annotation(
      Placement(transformation(origin = {-121, 113}, extent = {{-7, -7}, {7, 7}})));
    Modelica.Blocks.Sources.Ramp ramp(duration = 0, height = 3000, offset = 0, startTime = 10) annotation(
      Placement(transformation(origin = {-109, -5}, extent = {{-7, -7}, {7, 7}})));
    Modelica.Blocks.Sources.Ramp ramp1(duration = 1, height = 0.167, offset = 0, startTime = 0) annotation(
      Placement(transformation(origin = {-83, -61}, extent = {{-7, -7}, {7, 7}})));
    Buildings.Fluid.Sources.Boundary_ph bouOut(use_h_in = true, nPorts = 1, redeclare package Medium = Medium, use_X_in = false, use_Xi_in = false, use_C_in = false, use_p_in = false, p = 101325) annotation(
      Placement(transformation(origin = {42, -12}, extent = {{-6, -6}, {6, 6}}, rotation = -90)));
    Buildings.Fluid.Sources.Boundary_pT bouIn(nPorts = 1, p = 101325, T = 303.15, redeclare package Medium = Medium, use_Xi_in = false, use_C_in = false, use_p_in = false, use_T_in = false) annotation(
      Placement(transformation(origin = {26, -12}, extent = {{-6, -6}, {6, 6}}, rotation = 90)));
    Modelica.Fluid.Sensors.SpecificEnthalpyTwoPort specificEnthalpy(redeclare package Medium = Medium) annotation(
      Placement(transformation(origin = {26, 12}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  equation
    connect(pipe1.port_a, hex_hp.port_b) annotation(
      Line(points = {{-74, 36}, {-74, 26}}, color = {0, 127, 255}));
    connect(pipe1.port_b, pipe2.port_a) annotation(
      Line(points = {{-74, 56}, {-74, 96}, {-34, 96}}, color = {0, 127, 255}));
    connect(pipe4.port_b, hex_hp.port_a) annotation(
      Line(points = {{-74, -4}, {-74, 6}}, color = {0, 127, 255}));
    connect(pump_m_flow.port_b, pipe4.port_a) annotation(
      Line(points = {{-64, -40}, {-74, -40}, {-74, -24}}, color = {0, 127, 255}));
    connect(pipe1.port_b, val.port_a) annotation(
      Line(points = {{-74, 56}, {-74, 66}, {-58, 66}}, color = {0, 127, 255}));
    connect(val.port_b, pipe5.port_a) annotation(
      Line(points = {{-38, 66}, {-24, 66}, {-24, 52}}, color = {0, 127, 255}));
    connect(pipe5.port_b, pump_m_flow.port_a) annotation(
      Line(points = {{-24, 32}, {-24, -40}, {-44, -40}}, color = {0, 127, 255}));
    connect(const.y, val.y) annotation(
      Line(points = {{-56, 86}, {-48, 86}, {-48, 78}}, color = {0, 0, 127}));
    connect(outdoor.port, pipe2.heatPort) annotation(
      Line(points = {{-114, 113}, {-24, 113}, {-24, 106}}, color = {191, 0, 0}));
    connect(outdoor.port, pipe1.heatPort) annotation(
      Line(points = {{-114, 113}, {-90, 113}, {-90, 46}, {-84, 46}}, color = {191, 0, 0}));
    connect(outdoor.port, pipe4.heatPort) annotation(
      Line(points = {{-114, 113}, {-92, 113}, {-92, -14}, {-84, -14}}, color = {191, 0, 0}));
    connect(outdoor.port, pipe5.heatPort) annotation(
      Line(points = {{-114, 113}, {-6, 113}, {-6, 42}, {-14, 42}}, color = {191, 0, 0}));
    connect(outdoor.port, pipe3.heatPort) annotation(
      Line(points = {{-114, 113}, {46, 113}, {46, 40}, {36, 40}}, color = {191, 0, 0}));
    connect(ramp.y, hex_hp.freq) annotation(
      Line(points = {{-101, -5}, {-82, -5}, {-82, 6}}, color = {0, 0, 127}));
    connect(ramp1.y, pump_m_flow.m_flow_in) annotation(
      Line(points = {{-76, -60}, {-54, -60}, {-54, -52}}, color = {0, 0, 127}));
    connect(pipe3.port_b, specificEnthalpy.port_a) annotation(
      Line(points = {{26, 30}, {26, 22}}, color = {0, 127, 255}));
    connect(specificEnthalpy.port_b, bouIn.ports[1]) annotation(
      Line(points = {{26, 2}, {26, -6}}, color = {0, 127, 255}));
    connect(specificEnthalpy.h_out, bouOut.h_in) annotation(
      Line(points = {{38, 12}, {38, 11}, {44, 11}, {44, -5}}, color = {0, 0, 127}));
    connect(bouOut.ports[1], pump_m_flow.port_a) annotation(
      Line(points = {{42, -18}, {42, -40}, {-44, -40}}, color = {0, 127, 255}));
    connect(pipe2.port_b, BufferTank.port_a) annotation(
      Line(points = {{-14, 96}, {6, 96}, {6, 77}}, color = {0, 127, 255}));
    connect(BufferTank.port_b, pipe3.port_a) annotation(
      Line(points = {{36, 77}, {42, 77}, {42, 56}, {26, 56}, {26, 50}}, color = {0, 127, 255}));
    annotation(
      Diagram(coordinateSystem(extent = {{-140, 160}, {60, -80}})),
      experiment(StartTime = 0, StopTime = 10000, Tolerance = 1e-6, Interval = 20));
  end plantold1;

  model plantold2
    replaceable package Medium = Buildings.Media.Water;
    inner Modelica.Fluid.System system annotation(
      Placement(transformation(origin = {40, 142}, extent = {{-10, -10}, {10, 10}})));
    import SI = Modelica.Units.SI;
    // Parameters
    import Modelica.Constants.pi;
    parameter Modelica.Units.SI.Temperature Amb_T = 273.15 + 10 "Ambience Temperature [K]";
    // Initialization
    parameter Modelica.Units.SI.Temperature medium_initT = 273.15 + 30 "Water Initial Temperature [K]";
    parameter Modelica.Units.SI.Temperature BT_initT = 273.15 + 30 "Buffer Tank Init Temperature [K]";
    parameter Modelica.Units.SI.Pressure BT_initP = 101325 "Buffer Tank Init Pressure [Pa]";
    // Buffer Tank
    parameter Modelica.Units.SI.Volume BT_vol = 0.2 "Buffer Tank Volume [m3]";
    parameter Modelica.Units.SI.Length BT_height = 1 "Buffer Tank height [m]";
    // Pipe
    parameter Modelica.Units.SI.Length pip_len = 1 "Length of pipe [m]";
    parameter Modelica.Units.SI.SpecificHeatCapacity pip_c = 386 "Specific heat of pipe wall material [J/kg.K]";
    parameter Modelica.Units.SI.Density pip_rho = 8960 "Density of pipe wall material [kg/m3]";
    parameter Modelica.Units.SI.Length pip_dh = 0.025 "Inner diameter of pipe [m]";
    parameter Modelica.Units.SI.Length pip_thickness = 0.005 "Thickness of pipe wall [m]";
    parameter Modelica.Units.SI.Length pip_dIns = 0.005 "Thickness of pipe insulation [m]";
    //parameter Modelica.Units.SI.ThermalConductivity pip_kIns = 1 "Heat conductivity of pipe insulation [W/m.K]";
    parameter Modelica.Units.SI.ThermalConductivity pip_kIns = 398 "Heat conductivity of pipe insulation [W/m.K]";
    parameter Modelica.Units.SI.CoefficientOfHeatTransfer pip_air_htc = 10 "Air-Pipe heat transfer coeff. [W/(m2.K)]";
    // Pump
    parameter Real m_flow = 10/60*0.001*1000 "mass flow [kg/s]";
    //
    final parameter Modelica.Units.SI.MassFlowRate m_flow_nominal = 6 "Mass flow rate [kg/s]";
    final parameter Modelica.Units.SI.PressureDifference dp_nominal = 4500 "Design pressure drop";
    //===========
    // components
    Buildings.Fluid.Movers.FlowControlled_m_flow pump_m_flow(redeclare package Medium = Medium, m_flow_nominal = m_flow_nominal, use_inputFilter = false, massFlowRates = {0, 0.5, 1}*m_flow_nominal, inputType = Buildings.Fluid.Types.InputType.Continuous, energyDynamics = Modelica.Fluid.Types.Dynamics.SteadyState, nominalValuesDefineDefaultPressureCurve = true, allowFlowReversal = false, T_start = medium_initT) "Pump with m_flow input" annotation(
      Placement(transformation(origin = {-54, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
    MyComponents.HEX_HP hex_hp(InitialTemp = medium_initT, redeclare package Medium = Medium) "Heat Exchanger compresser frequency input" annotation(
      Placement(transformation(origin = {-74, 16}, extent = {{-12, -8}, {12, 14}}, rotation = 90)));
    Buildings.Fluid.FixedResistances.PlugFlowPipe pipe1(redeclare package Medium = Medium, length = pip_len, T_start_in = medium_initT, T_start_out = medium_initT, allowFlowReversal = false, cPip = pip_c, rhoPip = pip_rho, thickness = pip_thickness, dh = pip_dh, dIns = pip_dIns, kIns = pip_kIns, m_flow_nominal = m_flow_nominal) "Outside PlugFlowPipe" annotation(
      Placement(transformation(origin = {-74, 46}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
    Buildings.Fluid.Storage.Stratified BufferTank(redeclare package Medium = Medium, T_start = BT_initT, VTan = BT_vol, dIns = 0.01, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, hTan = BT_height, nSeg = 5, p_start = BT_initP) "comment" annotation(
      Placement(transformation(origin = {27, 81}, extent = {{15, -15}, {-15, 15}})));
    Buildings.Fluid.FixedResistances.PlugFlowPipe pipe3(redeclare package Medium = Medium, T_start_in = medium_initT, T_start_out = medium_initT, allowFlowReversal = false, cPip = pip_c, length = pip_len, rhoPip = pip_rho, thickness = pip_thickness, dh = pip_dh, dIns = pip_dIns, kIns = pip_kIns, m_flow_nominal = m_flow_nominal) "comment" annotation(
      Placement(transformation(origin = {26, 40}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Buildings.Fluid.FixedResistances.PlugFlowPipe pipe2(redeclare package Medium = Medium, T_start_in = medium_initT, T_start_out = medium_initT, allowFlowReversal = false, cPip = pip_c, length = pip_len, rhoPip = pip_rho, thickness = pip_thickness, dh = pip_dh, dIns = pip_dIns, kIns = pip_kIns, m_flow_nominal = m_flow_nominal) "Indoor PlugFlowPipe" annotation(
      Placement(transformation(origin = {-24, 96}, extent = {{-10, -10}, {10, 10}})));
    Buildings.Fluid.FixedResistances.PlugFlowPipe pipe4(redeclare package Medium = Medium, T_start_in = medium_initT, T_start_out = medium_initT, allowFlowReversal = false, cPip = pip_c, length = pip_len, rhoPip = pip_rho, thickness = pip_thickness, dh = pip_dh, dIns = pip_dIns, kIns = pip_kIns, m_flow_nominal = m_flow_nominal) "" annotation(
      Placement(transformation(origin = {-74, -14}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
    Buildings.Fluid.Actuators.Valves.TwoWayLinear val(redeclare package Medium = Medium, m_flow_nominal = m_flow_nominal, CvData = Buildings.Fluid.Types.CvTypes.OpPoint, dpValve_nominal(displayUnit = "kPa") = dp_nominal, use_inputFilter = false) "Valve model, linear opening characteristics" annotation(
      Placement(transformation(origin = {-48, 66}, extent = {{-10, -10}, {10, 10}})));
    Buildings.Fluid.FixedResistances.PlugFlowPipe pipe5(redeclare package Medium = Medium, T_start_in = medium_initT, T_start_out = medium_initT, allowFlowReversal = false, cPip = pip_c, length = pip_len, rhoPip = pip_rho, thickness = pip_thickness, dh = pip_dh, dIns = pip_dIns, kIns = pip_kIns, m_flow_nominal = m_flow_nominal) "" annotation(
      Placement(transformation(origin = {-24, 42}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Blocks.Sources.Constant const(k = 0) annotation(
      Placement(transformation(origin = {-62, 86}, extent = {{-6, -6}, {6, 6}})));
    Buildings.HeatTransfer.Sources.FixedTemperature outdoor(each T = Amb_T) "Boundary temperature" annotation(
      Placement(transformation(origin = {-140, 115}, extent = {{-7, -7}, {7, 7}})));
    Modelica.Blocks.Sources.Ramp ramp(duration = 0, height = 50, offset = 0, startTime = 10) annotation(
      Placement(transformation(origin = {-109, -5}, extent = {{-7, -7}, {7, 7}})));
    Modelica.Blocks.Sources.Ramp ramp1(duration = 1, height = m_flow, offset = 0, startTime = 0) annotation(
      Placement(transformation(origin = {-85, -61}, extent = {{-7, -7}, {7, 7}})));
    Buildings.Fluid.Storage.ExpansionVessel exp(redeclare package Medium = Medium, p_start = 101325) annotation(
      Placement(transformation(origin = {48, -20}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Thermal.HeatTransfer.Components.Convection convection annotation(
      Placement(transformation(origin = {-104, 115}, extent = {{7, -7}, {-7, 7}})));
    Modelica.Blocks.Sources.Constant Gc_pipe(k = pip_air_htc*pi*(pip_dh + pip_thickness)*pip_len) annotation(
      Placement(transformation(origin = {-119, 139}, extent = {{-7, -7}, {7, 7}})));
  equation
    connect(pipe1.port_a, hex_hp.port_b) annotation(
      Line(points = {{-74, 36}, {-74, 26}}, color = {0, 127, 255}));
    connect(pipe1.port_b, pipe2.port_a) annotation(
      Line(points = {{-74, 56}, {-74, 96}, {-34, 96}}, color = {0, 127, 255}));
    connect(pipe4.port_b, hex_hp.port_a) annotation(
      Line(points = {{-74, -4}, {-74, 6}}, color = {0, 127, 255}));
    connect(pump_m_flow.port_b, pipe4.port_a) annotation(
      Line(points = {{-64, -40}, {-74, -40}, {-74, -24}}, color = {0, 127, 255}));
    connect(pipe1.port_b, val.port_a) annotation(
      Line(points = {{-74, 56}, {-74, 66}, {-58, 66}}, color = {0, 127, 255}));
    connect(val.port_b, pipe5.port_a) annotation(
      Line(points = {{-38, 66}, {-24, 66}, {-24, 52}}, color = {0, 127, 255}));
    connect(pipe5.port_b, pump_m_flow.port_a) annotation(
      Line(points = {{-24, 32}, {-24, -40}, {-44, -40}}, color = {0, 127, 255}));
    connect(const.y, val.y) annotation(
      Line(points = {{-56, 86}, {-48, 86}, {-48, 78}}, color = {0, 0, 127}));
    connect(ramp.y, hex_hp.freq) annotation(
      Line(points = {{-101, -5}, {-82, -5}, {-82, 6}}, color = {0, 0, 127}));
    connect(ramp1.y, pump_m_flow.m_flow_in) annotation(
      Line(points = {{-77, -61}, {-77, -52}, {-54, -52}}, color = {0, 0, 127}));
    connect(BufferTank.port_b, pipe3.port_a) annotation(
      Line(points = {{12, 81}, {12, 56}, {26, 56}, {26, 50}}, color = {0, 127, 255}));
    connect(pipe2.port_b, BufferTank.port_a) annotation(
      Line(points = {{-14, 96}, {42, 96}, {42, 82}}, color = {0, 127, 255}));
    connect(pipe3.port_b, pump_m_flow.port_a) annotation(
      Line(points = {{26, 30}, {26, -40}, {-44, -40}}, color = {0, 127, 255}));
    connect(pipe3.port_b, exp.port_a) annotation(
      Line(points = {{26, 30}, {26, -20}, {38, -20}}, color = {0, 127, 255}));
    connect(outdoor.port, convection.fluid) annotation(
      Line(points = {{-133, 115}, {-111, 115}}, color = {191, 0, 0}));
    connect(Gc_pipe.y, convection.Gc) annotation(
      Line(points = {{-112, 140}, {-104, 140}, {-104, 122}}, color = {0, 0, 127}));
    connect(convection.solid, pipe2.heatPort) annotation(
      Line(points = {{-96, 116}, {-24, 116}, {-24, 106}}, color = {191, 0, 0}));
    connect(convection.solid, pipe5.heatPort) annotation(
      Line(points = {{-96, 116}, {-8, 116}, {-8, 42}, {-14, 42}}, color = {191, 0, 0}));
    connect(convection.solid, pipe3.heatPort) annotation(
      Line(points = {{-96, 116}, {52, 116}, {52, 40}, {36, 40}}, color = {191, 0, 0}));
    connect(convection.solid, pipe1.heatPort) annotation(
      Line(points = {{-96, 116}, {-88, 116}, {-88, 46}, {-84, 46}}, color = {191, 0, 0}));
    connect(convection.solid, pipe4.heatPort) annotation(
      Line(points = {{-96, 116}, {-92, 116}, {-92, -14}, {-84, -14}}, color = {191, 0, 0}));
    annotation(
      Diagram(coordinateSystem(extent = {{-140, 160}, {60, -80}})),
      experiment(StartTime = 0, StopTime = 100, Tolerance = 1e-06, Interval = 0.2));
  end plantold2;

  model plant
    replaceable package Medium = Buildings.Media.Water;
    inner Modelica.Fluid.System system annotation(
      Placement(transformation(origin = {40, 142}, extent = {{-10, -10}, {10, 10}})));
    import SI = Modelica.Units.SI;
    // Parameters
    import Modelica.Constants.pi;
    parameter Modelica.Units.SI.Temperature Amb_T = 273.15 + 5 "Ambience Temperature [K]";
    // Initialization
    parameter Modelica.Units.SI.Temperature medium_initT = 273.15 + 30 "Water Initial Temperature [K]";
    parameter Modelica.Units.SI.Temperature BT_initT = 273.15 + 30 "Buffer Tank Init Temperature [K]";
    parameter Modelica.Units.SI.Pressure BT_initP = 101325 "Buffer Tank Init Pressure [Pa]";
    // Buffer Tank
    parameter Modelica.Units.SI.Volume BT_vol = 0.2 "Buffer Tank Volume [m3]";
    parameter Modelica.Units.SI.Length BT_height = 1 "Buffer Tank height [m]";
    // Pipe
    parameter Modelica.Units.SI.Length pip_len = 1 "Length of pipe [m]";
    parameter Modelica.Units.SI.SpecificHeatCapacity pip_c = 386 "Specific heat of pipe wall material [J/kg.K]";
    parameter Modelica.Units.SI.Density pip_rho = 8960 "Density of pipe wall material [kg/m3]";
    parameter Modelica.Units.SI.Length pip_dh = 0.025 "Inner diameter of pipe [m]";
    parameter Modelica.Units.SI.Length pip_thickness = 0.005 "Thickness of pipe wall [m]";
    parameter Modelica.Units.SI.ThermalConductivity pip_kIns = 100 "Heat conductivity of pipe insulation [W/m.K]";
    parameter Modelica.Units.SI.CoefficientOfHeatTransfer pip_air_htc = 10 "Air-Pipe heat transfer coeff. [W/(m2.K)]";
    // Pump
    parameter Real m_flow = 10/60*0.001*1000 "mass flow [kg/s]";
    //
    final parameter Modelica.Units.SI.MassFlowRate m_flow_nominal = 6 "Mass flow rate [kg/s]";
    final parameter Modelica.Units.SI.PressureDifference dp_nominal = 4500 "Design pressure drop";
    //===========
    // components
    Buildings.Fluid.Movers.FlowControlled_m_flow pump_m_flow(redeclare package Medium = Medium, m_flow_nominal = m_flow_nominal, use_inputFilter = false, massFlowRates = {0, 0.5, 1}*m_flow_nominal, inputType = Buildings.Fluid.Types.InputType.Continuous, energyDynamics = Modelica.Fluid.Types.Dynamics.SteadyState, nominalValuesDefineDefaultPressureCurve = true, allowFlowReversal = false, T_start = medium_initT) "Pump with m_flow input" annotation(
      Placement(transformation(origin = {-54, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
    Buildings.Fluid.Storage.Stratified BufferTank(redeclare package Medium = Medium, T_start = BT_initT, VTan = BT_vol, dIns = 0.01, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, hTan = BT_height, nSeg = 2) "comment" annotation(
      Placement(transformation(origin = {27, 81}, extent = {{15, -15}, {-15, 15}})));
    Buildings.Fluid.FixedResistances.PlugFlowPipe pipe1(redeclare package Medium = Medium, length = pip_len, T_start_in = medium_initT, T_start_out = medium_initT, allowFlowReversal = false, cPip = pip_c, rhoPip = pip_rho, thickness = pip_thickness, dh = pip_dh, dIns = pip_thickness, kIns = pip_kIns, m_flow_nominal = m_flow_nominal, R = Modelica.Math.log((pip_dh + pip_thickness)/pip_dh)/(2*Modelica.Constants.pi*pip_kIns) + 1/pip_air_htc/Modelica.Constants.pi/(pip_dh + pip_thickness)) "Outside PlugFlowPipe" annotation(
      Placement(transformation(origin = {-74, 46}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
    Buildings.Fluid.FixedResistances.PlugFlowPipe pipe2(redeclare package Medium = Medium, T_start_in = medium_initT, T_start_out = medium_initT, allowFlowReversal = false, cPip = pip_c, length = pip_len, rhoPip = pip_rho, thickness = pip_thickness, dh = pip_dh, dIns = pip_thickness, kIns = pip_kIns, m_flow_nominal = m_flow_nominal, R = Modelica.Math.log((pip_dh + pip_thickness)/pip_dh)/(2*Modelica.Constants.pi*pip_kIns) + 1/pip_air_htc/Modelica.Constants.pi/(pip_dh + pip_thickness)) "Indoor PlugFlowPipe" annotation(
      Placement(transformation(origin = {-24, 96}, extent = {{-10, -10}, {10, 10}})));
    Buildings.Fluid.FixedResistances.PlugFlowPipe pipe3(redeclare package Medium = Medium, T_start_in = medium_initT, T_start_out = medium_initT, allowFlowReversal = false, cPip = pip_c, length = pip_len, rhoPip = pip_rho, thickness = pip_thickness, dh = pip_dh, dIns = pip_thickness, kIns = pip_kIns, m_flow_nominal = m_flow_nominal, R = Modelica.Math.log((pip_dh + pip_thickness)/pip_dh)/(2*Modelica.Constants.pi*pip_kIns) + 1/pip_air_htc/Modelica.Constants.pi/(pip_dh + pip_thickness)) "comment" annotation(
      Placement(transformation(origin = {26, 40}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Buildings.Fluid.FixedResistances.PlugFlowPipe pipe4(redeclare package Medium = Medium, T_start_in = medium_initT, T_start_out = medium_initT, allowFlowReversal = false, cPip = pip_c, length = pip_len, rhoPip = pip_rho, thickness = pip_thickness, dh = pip_dh, dIns = pip_thickness, kIns = pip_kIns, m_flow_nominal = m_flow_nominal, R = Modelica.Math.log((pip_dh + pip_thickness)/pip_dh)/(2*Modelica.Constants.pi*pip_kIns) + 1/pip_air_htc/Modelica.Constants.pi/(pip_dh + pip_thickness)) "" annotation(
      Placement(transformation(origin = {-74, -14}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
    Buildings.Fluid.FixedResistances.PlugFlowPipe pipe5(redeclare package Medium = Medium, T_start_in = medium_initT, T_start_out = medium_initT, allowFlowReversal = false, cPip = pip_c, length = pip_len, rhoPip = pip_rho, thickness = pip_thickness, dh = pip_dh, dIns = pip_thickness, kIns = pip_kIns, m_flow_nominal = m_flow_nominal, R = Modelica.Math.log((pip_dh + pip_thickness)/pip_dh)/(2*Modelica.Constants.pi*pip_kIns) + 1/pip_air_htc/Modelica.Constants.pi/(pip_dh + pip_thickness)) "" annotation(
      Placement(transformation(origin = {-24, 42}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Buildings.Fluid.Actuators.Valves.TwoWayLinear val(redeclare package Medium = Medium, m_flow_nominal = m_flow_nominal, CvData = Buildings.Fluid.Types.CvTypes.OpPoint, dpValve_nominal(displayUnit = "kPa") = dp_nominal, use_inputFilter = false) "Valve model, linear opening characteristics" annotation(
      Placement(transformation(origin = {-48, 66}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Sources.Constant const(k = 0) annotation(
      Placement(transformation(origin = {-62, 86}, extent = {{-6, -6}, {6, 6}})));
    Buildings.HeatTransfer.Sources.FixedTemperature outdoor(each T = Amb_T) "Boundary temperature" annotation(
      Placement(transformation(origin = {-118, 117}, extent = {{-7, -7}, {7, 7}})));
    Modelica.Blocks.Sources.Ramp ramp(duration = 0, height = 80, offset = 0, startTime = 3600) annotation(
      Placement(transformation(origin = {-109, -5}, extent = {{-7, -7}, {7, 7}})));
    Modelica.Blocks.Sources.Ramp ramp1(duration = 0, height = m_flow, offset = 0, startTime = 0) annotation(
      Placement(transformation(origin = {-85, -61}, extent = {{-7, -7}, {7, 7}})));
    Buildings.Fluid.Storage.ExpansionVessel exp(redeclare package Medium = Medium, p_start = 101325, V_start = 0.001, T_start = medium_initT) annotation(
      Placement(transformation(origin = {26, -54}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
    Modelica.Fluid.Sensors.TemperatureTwoPort temperature_pipe2(redeclare package Medium = Medium) annotation(
      Placement(transformation(origin = {18, 110}, extent = {{-6, -6}, {6, 6}})));
    Modelica.Fluid.Sensors.TemperatureTwoPort temperature_pipe1(redeclare package Medium = Medium) annotation(
      Placement(transformation(origin = {-50, 106}, extent = {{-6, -6}, {6, 6}})));
    MyComponents.roomCycleWithoutLossnay room1(q_flow_nominal = 2500) annotation(
      Placement(transformation(origin = {74, 76}, extent = {{-10, -10}, {10, 10}})));
    MyComponents.roomCycleWithoutLossnay room2(room_w = 3.2, q_flow_nominal = 3750) annotation(
      Placement(transformation(origin = {74, 52}, extent = {{-10, -10}, {10, 10}})));
    MyComponents.roomCycleWithoutLossnay room3(room_w = 4.8, q_flow_nominal = 5000) annotation(
      Placement(transformation(origin = {74, 28}, extent = {{-10, -10}, {10, 10}})));
    MyComponents.roomCycleWithoutLossnay room4(room_w = 6.4, q_flow_nominal = 6250) annotation(
      Placement(transformation(origin = {74, 4}, extent = {{-10, -10}, {10, 10}})));
    MyComponents.roomCycleWithoutLossnay room5(room_w = 8.0, q_flow_nominal = 7500) annotation(
      Placement(transformation(origin = {74, -20}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Thermal.HeatTransfer.Celsius.TemperatureSensor temperatureSensor1 annotation(
      Placement(transformation(origin = {96, 100}, extent = {{-10, -10}, {10, 10}})));
    System.MyComponents.HEX_HP hex_hp(InitialTemp = medium_initT) annotation(
      Placement(visible = true, transformation(origin = {-77.6112, 18.3373}, extent = {{-14.6707, 15.2349}, {7.33533, 37.2409}}, rotation = 90)));
  equation
    connect(pump_m_flow.port_b, pipe4.port_a) annotation(
      Line(points = {{-64, -40}, {-74, -40}, {-74, -24}}, color = {0, 127, 255}));
    connect(pipe1.port_b, val.port_a) annotation(
      Line(points = {{-74, 56}, {-74, 66}, {-58, 66}}, color = {0, 127, 255}));
    connect(val.port_b, pipe5.port_a) annotation(
      Line(points = {{-38, 66}, {-24, 66}, {-24, 52}}, color = {0, 127, 255}));
    connect(pipe5.port_b, pump_m_flow.port_a) annotation(
      Line(points = {{-24, 32}, {-24, -40}, {-44, -40}}, color = {0, 127, 255}));
    connect(const.y, val.y) annotation(
      Line(points = {{-56, 86}, {-48, 86}, {-48, 78}}, color = {0, 0, 127}));
    connect(ramp1.y, pump_m_flow.m_flow_in) annotation(
      Line(points = {{-77, -61}, {-77, -52}, {-54, -52}}, color = {0, 0, 127}));
    connect(BufferTank.port_b, pipe3.port_a) annotation(
      Line(points = {{12, 81}, {12, 56}, {26, 56}, {26, 50}}, color = {0, 127, 255}));
    connect(pipe3.port_b, pump_m_flow.port_a) annotation(
      Line(points = {{26, 30}, {26, -40}, {-44, -40}}, color = {0, 127, 255}));
    connect(pipe3.port_b, exp.port_a) annotation(
      Line(points = {{26, 30}, {26, -44}}, color = {0, 127, 255}));
    connect(outdoor.port, pipe2.heatPort) annotation(
      Line(points = {{-110, 118}, {-24, 118}, {-24, 106}}, color = {191, 0, 0}));
    connect(outdoor.port, pipe3.heatPort) annotation(
      Line(points = {{-110, 118}, {50, 118}, {50, 40}, {36, 40}}, color = {191, 0, 0}));
    connect(outdoor.port, pipe5.heatPort) annotation(
      Line(points = {{-110, 118}, {-4, 118}, {-4, 42}, {-14, 42}}, color = {191, 0, 0}));
    connect(outdoor.port, pipe1.heatPort) annotation(
      Line(points = {{-110, 118}, {-92, 118}, {-92, 46}, {-84, 46}}, color = {191, 0, 0}));
    connect(outdoor.port, pipe4.heatPort) annotation(
      Line(points = {{-110, 118}, {-92, 118}, {-92, -14}, {-84, -14}}, color = {191, 0, 0}));
    connect(pipe2.port_b, temperature_pipe2.port_a) annotation(
      Line(points = {{-14, 96}, {4, 96}, {4, 110}, {12, 110}}, color = {0, 127, 255}));
    connect(temperature_pipe2.port_b, BufferTank.port_a) annotation(
      Line(points = {{24, 110}, {48, 110}, {48, 82}, {42, 82}}, color = {0, 127, 255}));
    connect(temperature_pipe1.port_b, pipe2.port_a) annotation(
      Line(points = {{-44, 106}, {-38, 106}, {-38, 96}, {-34, 96}}, color = {0, 127, 255}));
    connect(pipe1.port_b, temperature_pipe1.port_a) annotation(
      Line(points = {{-74, 56}, {-74, 106}, {-56, 106}}, color = {0, 127, 255}));
    connect(BufferTank.fluPorVol[1], room1.port_a) annotation(
      Line(points = {{30, 82}, {64, 82}, {64, 81}}, color = {0, 127, 255}));
    connect(BufferTank.fluPorVol[2], room1.port_b) annotation(
      Line(points = {{30, 82}, {32, 82}, {32, 71}, {64, 71}}, color = {0, 127, 255}));
    connect(BufferTank.fluPorVol[1], room2.port_a) annotation(
      Line(points = {{30, 82}, {32, 82}, {32, 58}, {64, 58}}, color = {0, 127, 255}));
    connect(BufferTank.fluPorVol[2], room2.port_b) annotation(
      Line(points = {{30, 82}, {32, 82}, {32, 48}, {64, 48}}, color = {0, 127, 255}));
    connect(BufferTank.fluPorVol[1], room3.port_a) annotation(
      Line(points = {{30, 82}, {32, 82}, {32, 34}, {64, 34}}, color = {0, 127, 255}));
    connect(BufferTank.fluPorVol[2], room3.port_b) annotation(
      Line(points = {{30, 82}, {32, 82}, {32, 24}, {64, 24}}, color = {0, 127, 255}));
    connect(BufferTank.fluPorVol[1], room4.port_a) annotation(
      Line(points = {{30, 82}, {32, 82}, {32, 10}, {64, 10}}, color = {0, 127, 255}));
    connect(BufferTank.fluPorVol[2], room4.port_b) annotation(
      Line(points = {{30, 82}, {32, 82}, {32, -2}, {64, -2}}, color = {0, 127, 255}));
    connect(BufferTank.fluPorVol[1], room5.port_a) annotation(
      Line(points = {{30, 82}, {32, 82}, {32, -14}, {64, -14}}, color = {0, 127, 255}));
    connect(BufferTank.fluPorVol[2], room5.port_b) annotation(
      Line(points = {{30, 82}, {32, 82}, {32, -24}, {64, -24}}, color = {0, 127, 255}));
    connect(BufferTank.heaPorVol[1], temperatureSensor1.port) annotation(
      Line(points = {{28, 82}, {86, 82}, {86, 100}}, color = {191, 0, 0}));
  connect(pipe1.port_a, hex_hp.port_b) annotation(
      Line(points = {{-74, 36}, {-74, 32}, {-73, 32}, {-73, 26}}, color = {0, 127, 255}));
  connect(pipe4.port_b, hex_hp.port_a) annotation(
      Line(points = {{-74, -4}, {-74, 2}, {-73, 2}, {-73, 11}}, color = {0, 127, 255}));
  connect(ramp.y, hex_hp.freq) annotation(
      Line(points = {{-101, -5}, {-101, 1.5}, {-83, 1.5}, {-83, 8}}, color = {0, 0, 127}));
    annotation(
      Diagram(coordinateSystem(extent = {{-140, 160}, {80, -80}})),
      experiment(StartTime = 0, StopTime = 14400, Tolerance = 1e-06, Interval = 28.8));
  end plant;
  annotation(
    uses(Modelica(version = "4.0.0"), Buildings(version = "10.0.0")));
end System;