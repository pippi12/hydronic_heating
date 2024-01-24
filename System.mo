package System
  package MyComponents
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
      Buildings.Fluid.Movers.FlowControlled_m_flow pump_m_flow_SA(redeclare package Medium = MediumAir, T_start = mediumRoomAir_initT, allowFlowReversal = true, energyDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial, inputType = Buildings.Fluid.Types.InputType.Continuous, m_flow_nominal = m_flow_nominal, massFlowRates = {0, 0.5, 1} * m_flow_nominal, nominalValuesDefineDefaultPressureCurve = true, use_inputFilter = false) annotation(
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
      Buildings.Fluid.Movers.FlowControlled_m_flow pump_m_flow_RA(redeclare package Medium = MediumAir, T_start = mediumRoomAir_initT, allowFlowReversal = true, energyDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial, inputType = Buildings.Fluid.Types.InputType.Continuous, m_flow_nominal = m_flow_nominal, massFlowRates = {0, 0.5, 1} * m_flow_nominal, nominalValuesDefineDefaultPressureCurve = true, use_inputFilter = false) annotation(
        Placement(transformation(origin = {132, -16}, extent = {{-10, -10}, {10, 10}})));
      Buildings.Fluid.Movers.FlowControlled_m_flow pump_m_flow_SA(redeclare package Medium = MediumAir, T_start = mediumRoomAir_initT, allowFlowReversal = true, energyDynamics = Modelica.Fluid.Types.Dynamics.DynamicFreeInitial, inputType = Buildings.Fluid.Types.InputType.Continuous, m_flow_nominal = m_flow_nominal, massFlowRates = {0, 0.5, 1} * m_flow_nominal, nominalValuesDefineDefaultPressureCurve = true, use_inputFilter = false) annotation(
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

    model HeatPump
      import SI = Modelica.Units.SI;
      replaceable package Medium = Buildings.Media.Water "Medium in the pipe";
      inner Modelica.Fluid.System system annotation(
        Placement(visible = true, transformation(origin = {-128, 88}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      // parameters
      parameter Real K = 150 "1st delay K";
      parameter Modelica.Units.SI.Time T = 500 "1st delay T";
      parameter Modelica.Units.SI.Time L = 100 "wasted time L";
      parameter Modelica.Units.SI.Temperature InitialTemp = 273.15 + 30 "Initial Temperature";
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
        Icon(graphics = {Rectangle(origin = {-20, 50}, fillColor = {77, 77, 77}, fillPattern = FillPattern.Solid, extent = {{-100, 90}, {100, -90}}), Text(origin = {-20, -60}, lineColor = {0, 0, 255}, extent = {{-160, 20}, {160, -20}}, textString = "%name"), Rectangle(origin = {-70, 0}, lineColor = {0, 0, 255}, fillColor = {0, 0, 255}, fillPattern = FillPattern.Solid, extent = {{-50, 4}, {50, -4}}), Rectangle(origin = {30, 0}, lineColor = {255, 0, 0}, fillColor = {255, 0, 0}, fillPattern = FillPattern.Solid, extent = {{-50, 4}, {50, -4}}), Text(origin = {-160, 155}, extent = {{-26, 27}, {26, -27}}, textString = "f"), Line(origin = {-130, 118}, points = {{-10, 0}, {10, 0}, {10, 0}}, thickness = 0.5)}, coordinateSystem(extent = {{-200, 180}, {100, -80}})),
        Diagram(coordinateSystem(extent = {{-140, -140}, {140, 140}})));
    end HeatPump;

    model CompressorController
      //parameters
      parameter Real dT = 1.5 "Temperature difference from OFF to ON";
      parameter Real lower = 20 "Lower limit frequency";
      parameter Real upper = 80 "Upper limit frequency";
      Modelica.Blocks.Interfaces.RealInput T_set annotation(
        Placement(visible = true, transformation(origin = {-102, 42}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-100, 60}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealInput T_BT annotation(
        Placement(visible = true, transformation(origin = {-104, -28}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-100, -60}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      Modelica.Blocks.Nonlinear.Limiter limiterMax(uMax = upper, uMin = 0) annotation(
        Placement(visible = true, transformation(origin = {34, 42}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Math.Feedback feedback1 annotation(
        Placement(visible = true, transformation(origin = {-40, 42}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Continuous.PI pi(T = 500, initType = Modelica.Blocks.Types.Init.InitialState, k = 5, x_start = 0, y_start = 0) annotation(
        Placement(visible = true, transformation(origin = {-2, 42}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Logical.Switch switch1 annotation(
        Placement(visible = true, transformation(origin = {108, 42}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression off_freq(y = 0) annotation(
        Placement(visible = true, transformation(origin = {74, 24}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.StateGraph.InitialStep initialStep(nIn = 1, nOut = 1) annotation(
        Placement(visible = true, transformation(extent = {{-32, 132}, {-12, 152}}, rotation = 0)));
      Modelica.StateGraph.StepWithSignal step_on_off(nIn = 1, nOut = 1) annotation(
        Placement(visible = true, transformation(extent = {{24, 132}, {44, 152}}, rotation = 0)));
      Modelica.StateGraph.TransitionWithSignal transition2 annotation(
        Placement(visible = true, transformation(extent = {{54, 132}, {74, 152}}, rotation = 0)));
      inner Modelica.StateGraph.StateGraphRoot stateGraphRoot annotation(
        Placement(visible = true, transformation(extent = {{-70, 136}, {-50, 156}}, rotation = 0)));
      Modelica.StateGraph.TransitionWithSignal transition1 annotation(
        Placement(visible = true, transformation(extent = {{-6, 132}, {14, 152}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealOutput y annotation(
        Placement(visible = true, transformation(origin = {148, 42}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {104, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Logical.GreaterThreshold greaterThreshold(threshold = dT) annotation(
        Placement(visible = true, transformation(origin = {-22, 74}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
      Modelica.Blocks.Logical.Not not1 annotation(
        Placement(visible = true, transformation(origin = {20, 86}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Logical.And and1 annotation(
        Placement(visible = true, transformation(origin = {64, 112}, extent = {{-10, 10}, {10, -10}}, rotation = 90)));
      Modelica.Blocks.Logical.LessEqualThreshold lessEqualThreshold(threshold = lower) annotation(
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
        Icon(graphics = {Text(origin = {-48, -60}, extent = {{-32, 12}, {32, -12}}, textString = "T_BT"), Text(origin = {-42, 60}, extent = {{-34, 12}, {34, -12}}, textString = "T_set"), Rectangle(extent = {{-100, 100}, {100, -100}}), Text(origin = {72, -1}, extent = {{-18, 17}, {18, -17}}, textString = "f"), Text(origin = {0, 120}, lineColor = {0, 0, 255}, extent = {{-100, 20}, {100, -20}}, textString = "%name")}),
        Diagram(coordinateSystem(extent = {{-120, 160}, {160, -60}})));
    end CompressorController;

    model PumpController ""
      //parameters
      parameter Real m_flow = 5 / 60 * 0.001 * 1000 "mass flow rate (kg/s)";
      parameter Real lower = -1 "Lower limit hysteresis (Tset-T)";
      parameter Real upper = 1 "Upper limit hysteresis (Tset-T)";
      Modelica.Blocks.Interfaces.RealInput T_set annotation(
        Placement(visible = true, transformation(origin = {-86, 42}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-100, 60}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealInput T_Air annotation(
        Placement(visible = true, transformation(origin = {-86, 12}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-100, -60}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      Modelica.Blocks.Math.Feedback feedback1 annotation(
        Placement(visible = true, transformation(origin = {-40, 42}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Logical.Switch switch annotation(
        Placement(visible = true, transformation(origin = {70, 42}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression off(y = 0) annotation(
        Placement(visible = true, transformation(origin = {36, 34}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealOutput y annotation(
        Placement(visible = true, transformation(origin = {102, 42}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {104, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Logical.Hysteresis hysteresis(pre_y_start = false, uHigh = 1, uLow = -1) annotation(
        Placement(visible = true, transformation(origin = {-6, 42}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression on(y = m_flow) annotation(
        Placement(visible = true, transformation(origin = {36, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(T_Air, feedback1.u2) annotation(
        Line(points = {{-86, 12}, {-40, 12}, {-40, 34}}, color = {0, 0, 127}));
      connect(T_set, feedback1.u1) annotation(
        Line(points = {{-86, 42}, {-48, 42}}, color = {0, 0, 127}));
      connect(switch.y, y) annotation(
        Line(points = {{81, 42}, {102, 42}}, color = {0, 0, 127}));
      connect(off.y, switch.u3) annotation(
        Line(points = {{47, 34}, {58, 34}}, color = {0, 0, 127}));
      connect(feedback1.y, hysteresis.u) annotation(
        Line(points = {{-30, 42}, {-18, 42}}, color = {0, 0, 127}));
      connect(hysteresis.y, switch.u2) annotation(
        Line(points = {{6, 42}, {58, 42}}, color = {255, 0, 255}));
      connect(on.y, switch.u1) annotation(
        Line(points = {{47, 50}, {57, 50}}, color = {0, 0, 127}));
      annotation(
        Icon(graphics = {Text(origin = {-48, -60}, extent = {{-32, 12}, {32, -12}}, textString = "T_Air"), Text(origin = {-42, 60}, extent = {{-34, 12}, {34, -12}}, textString = "T_set"), Rectangle(extent = {{-100, 100}, {100, -100}}), Text(origin = {58, -1}, extent = {{-32, 17}, {32, -17}}, textString = "m_flow"), Text(origin = {0, 120}, lineColor = {0, 0, 255}, extent = {{-100, 18}, {100, -18}}, textString = "%name")}),
        Diagram(coordinateSystem(extent = {{-120, 60}, {120, -20}})));
    end PumpController;

    model RoomWithoutLossnay
      replaceable package Medium = Buildings.Media.Air;
      inner Modelica.Fluid.System system annotation(
        Placement(transformation(origin = {-154, 102}, extent = {{-10, -10}, {10, 10}})));
      // Parameters
      parameter Modelica.Units.SI.Temperature amb_T = 273.15 + 5 "Ambience initial temperature";
      // Initialization
      parameter Modelica.Units.SI.Temperature mediumRoomAir_initT = 273.15 + 20 "Room air initial temperature";
      parameter Modelica.Units.SI.Pressure mediumRoomAir_initP = 101325 "Room air initial Pressure";
      parameter Modelica.Units.SI.Temperature rw_initT = 273.15 + 15 "Room wall initial temperature";
      // Room
      parameter Modelica.Units.SI.Length room_w = 1.6 "Width of room";
      parameter Modelica.Units.SI.Length room_d = 4.8 "Depth of room";
      parameter Modelica.Units.SI.Length room_h = 2.4 "Height of room";
      parameter Modelica.Units.SI.Area room_area = 2 * (room_w * room_d + room_d * room_h + room_h * room_w) "Area of room";
      parameter Modelica.Units.SI.Length rw_thickness = 0.1 "Thickness of room wall";
      parameter Modelica.Units.SI.Density rw_rho = 144 "Density of room wall";
      parameter Modelica.Units.SI.SpecificHeatCapacity rw_c = 1168 "Specific heat capacity of room wall";
      parameter Modelica.Units.SI.ThermalConductivity rw_k = 0.3 "Heat conductivity of room wall";
      parameter Modelica.Units.SI.Volume rw_vol = room_area * rw_thickness;
      parameter Modelica.Units.SI.CoefficientOfHeatTransfer rw_air_htc = 5 "Air-Room wall heat transfer coeff.";
      final parameter Modelica.Units.SI.MassFlowRate m_flow_nominal = 0.1 "Mass flow rate";
      //final parameter Modelica.Units.SI.PressureDifference dp_nominal = 4500 "Design pressure drop";
      // Components
      //
      Buildings.HeatTransfer.Sources.FixedTemperature ambT(T = amb_T) annotation(
        Placement(visible = true, transformation(origin = {-160, 34}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
      Buildings.Fluid.MixingVolumes.MixingVolume roomAir(redeclare package Medium = Medium, T_start = mediumRoomAir_initT, V = room_w * room_d * room_h, allowFlowReversal = true, m_flow_nominal = 0.01, massDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, nPorts = 0, p_start(displayUnit = "Pa") = mediumRoomAir_initP, use_C_flow = false) annotation(
        Placement(visible = true, transformation(origin = {22, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.HeatTransfer.Components.HeatCapacitor heatCapacitor(C = rw_vol * rw_rho * rw_c, T(fixed = true, start = rw_initT)) annotation(
        Placement(visible = true, transformation(origin = {-68, 54}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor_wall1(G = rw_k * room_area * 2 / rw_thickness) annotation(
        Placement(visible = true, transformation(origin = {-92, 34}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor_wall2(G = rw_k * room_area * 2 / rw_thickness) annotation(
        Placement(visible = true, transformation(origin = {-44, 34}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      Buildings.HeatTransfer.Convection.Exterior con_ext(A = room_area, azi = 3.14 / 2, conMod = Buildings.HeatTransfer.Types.ExteriorConvection.Fixed, hFixed = rw_air_htc, til = 0) annotation(
        Placement(visible = true, transformation(origin = {-132, 34}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      Buildings.HeatTransfer.Convection.Interior con_int(A = room_area, conMod = Buildings.HeatTransfer.Types.InteriorConvection.Fixed, hFixed = rw_air_htc, til = 0) annotation(
        Placement(visible = true, transformation(origin = {-6, 34}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant const_v(k = 0) annotation(
        Placement(visible = true, transformation(origin = {-94, 84}, extent = {{-7, 7}, {7, -7}}, rotation = -180)));
      Modelica.Blocks.Sources.Constant const_d(k = 0) annotation(
        Placement(visible = true, transformation(origin = {-94, 64}, extent = {{-7, 7}, {7, -7}}, rotation = -180)));
      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a HeatPortCon annotation(
        Placement(visible = true, transformation(origin = {-182, 4}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-2, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a HeatPortRad annotation(
        Placement(visible = true, transformation(origin = {-182, -22}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-2, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.HeatTransfer.Celsius.TemperatureSensor temperatureSensor annotation(
        Placement(visible = true, transformation(origin = {58, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealOutput T annotation(
        Placement(visible = true, transformation(origin = {98, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {150, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(thermalConductor_wall1.port_b, heatCapacitor.port) annotation(
        Line(points = {{-82, 34}, {-74, 34}, {-74, 44}, {-68, 44}}, color = {191, 0, 0}));
      connect(heatCapacitor.port, thermalConductor_wall2.port_b) annotation(
        Line(points = {{-68, 44}, {-62, 44}, {-62, 34}, {-54, 34}}, color = {191, 0, 0}));
      connect(con_ext.solid, thermalConductor_wall1.port_a) annotation(
        Line(points = {{-122, 34}, {-102, 34}}, color = {191, 0, 0}));
      connect(con_ext.fluid, ambT.port) annotation(
        Line(points = {{-142, 34}, {-152, 34}}, color = {191, 0, 0}));
      connect(thermalConductor_wall2.port_a, con_int.solid) annotation(
        Line(points = {{-34, 34}, {-16, 34}}, color = {191, 0, 0}));
      connect(con_int.fluid, roomAir.heatPort) annotation(
        Line(points = {{4, 34}, {8, 34}, {8, 2}, {12, 2}}, color = {191, 0, 0}));
      connect(const_v.y, con_ext.v) annotation(
        Line(points = {{-101.7, 84}, {-111.7, 84}, {-111.7, 44}, {-119.7, 44}}, color = {0, 0, 127}));
      connect(const_d.y, con_ext.dir) annotation(
        Line(points = {{-101.7, 64}, {-109.7, 64}, {-109.7, 40}, {-119.7, 40}}, color = {0, 0, 127}));
      connect(HeatPortCon, roomAir.heatPort) annotation(
        Line(points = {{-182, 4}, {12, 4}, {12, 2}}, color = {191, 0, 0}));
      connect(HeatPortRad, roomAir.heatPort) annotation(
        Line(points = {{-182, -22}, {-10, -22}, {-10, 2}, {12, 2}}, color = {191, 0, 0}));
      connect(roomAir.heatPort, temperatureSensor.port) annotation(
        Line(points = {{12, 2}, {4, 2}, {4, -20}, {48, -20}, {48, 2}}, color = {191, 0, 0}));
      connect(T, temperatureSensor.T) annotation(
        Line(points = {{98, 2}, {68, 2}}, color = {0, 0, 127}));
/*
      for i in 1:nPorts loop
        connect(ports[i],roomAir.ports[i])
        annotation (Line(
          points={{-260,-60},{-218,-60},{-218,-206},{52,-206},{52,-141.9}},
          color={0,127,255},
          smooth=Smooth.None));
      end for;
      */
      annotation(
        Icon(graphics = {Rectangle(origin = {0, 110}, fillColor = {103, 103, 103}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, lineThickness = 0, extent = {{-100, 10}, {100, -10}}), Rectangle(lineColor = {0, 85, 255}, fillColor = {0, 170, 255}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-100, 100}, {100, -100}}), Text(origin = {2, 133}, lineColor = {0, 0, 255}, extent = {{-104, 21}, {104, -21}}, textString = "%name"), Rectangle(origin = {0, -110}, fillColor = {103, 103, 103}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, lineThickness = 0, extent = {{-100, 10}, {100, -10}}), Rectangle(origin = {-110, 0}, fillColor = {103, 103, 103}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, lineThickness = 0, extent = {{-10, 120}, {10, -120}}), Rectangle(origin = {110, 0}, fillColor = {103, 103, 103}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, lineThickness = 0, extent = {{-10, 120}, {10, -120}}), Line(origin = {116.911, -0.117096}, points = {{-23, 0}, {23, 0}, {23, 0}})}, coordinateSystem(extent = {{-120, 160}, {160, -120}})),
        Diagram(coordinateSystem(extent = {{-200, 120}, {120, -40}})));
    end RoomWithoutLossnay;

    model RoomWithLossnay
      replaceable package Medium = Buildings.Media.Air;
      inner Modelica.Fluid.System system annotation(
        Placement(visible = true, transformation(origin = {-180, 76}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      // Parameters
      parameter Modelica.Units.SI.Pressure amb_P = 101325 "Ambience pressure";
      // Initialization
      parameter Modelica.Units.SI.Temperature amb_T = 273.15 + 5 "Ambience initial temperature";
      parameter Modelica.Units.SI.Temperature mediumRoomAir_initT = 273.15 + 20 "Room air initial temperature";
      parameter Modelica.Units.SI.Pressure mediumRoomAir_initP = 101325 "Room air initial Pressure";
      parameter Modelica.Units.SI.Temperature rw_initT = 273.15 + 15 "Room wall initial temperature";
      // Room
      parameter Modelica.Units.SI.Length room_w = 1.6 "Width of room";
      parameter Modelica.Units.SI.Length room_d = 4.8 "Depth of room";
      parameter Modelica.Units.SI.Length room_h = 2.4 "Height of room";
      parameter Modelica.Units.SI.Area room_area = 2 * (room_w * room_d + room_d * room_h + room_h * room_w) "Area of room";
      parameter Modelica.Units.SI.Length rw_thickness = 0.1 "Thickness of room wall";
      parameter Modelica.Units.SI.Density rw_rho = 144 "Density of room wall";
      parameter Modelica.Units.SI.SpecificHeatCapacity rw_c = 1168 "Specific heat capacity of room wall";
      parameter Modelica.Units.SI.ThermalConductivity rw_k = 0.3 "Heat conductivity of room wall";
      parameter Modelica.Units.SI.Volume rw_vol = room_area * rw_thickness;
      parameter Modelica.Units.SI.CoefficientOfHeatTransfer rw_air_htc = 5 "Air-Room wall heat transfer coeff.";
      final parameter Modelica.Units.SI.MassFlowRate m_flow_nominal = 0.1 "Mass flow rate";
      final parameter Modelica.Units.SI.MassFlowRate m_flow_nominal_lossnay = 2 * (room_w * room_d * room_h) * 1.2 * 0.37 / 3600 "Mass flow rate (Lossnay)";
      //final parameter Modelica.Units.SI.PressureDifference dp_nominal = 4500 "Design pressure drop";
      // Components
      //
      Buildings.HeatTransfer.Sources.FixedTemperature ambT(T = amb_T) annotation(
        Placement(visible = true, transformation(origin = {-160, 34}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
      Buildings.Fluid.MixingVolumes.MixingVolume roomAir(redeclare package Medium = Medium, T_start = mediumRoomAir_initT, V = room_w * room_d * room_h, allowFlowReversal = true, m_flow_nominal = 0.01, massDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, nPorts = 2, p_start(displayUnit = "Pa") = mediumRoomAir_initP, use_C_flow = false) annotation(
        Placement(visible = true, transformation(origin = {22, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.HeatTransfer.Components.HeatCapacitor heatCapacitor(C = rw_vol * rw_rho * rw_c, T(fixed = true, start = rw_initT)) annotation(
        Placement(visible = true, transformation(origin = {-68, 54}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor_wall1(G = rw_k * room_area * 2 / rw_thickness) annotation(
        Placement(visible = true, transformation(origin = {-92, 34}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor_wall2(G = rw_k * room_area * 2 / rw_thickness) annotation(
        Placement(visible = true, transformation(origin = {-44, 34}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      Buildings.HeatTransfer.Convection.Exterior con_ext(A = room_area, azi = 3.14 / 2, conMod = Buildings.HeatTransfer.Types.ExteriorConvection.Fixed, hFixed = rw_air_htc, til = 0) annotation(
        Placement(visible = true, transformation(origin = {-132, 34}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      Buildings.HeatTransfer.Convection.Interior con_int(A = room_area, conMod = Buildings.HeatTransfer.Types.InteriorConvection.Fixed, hFixed = rw_air_htc, til = 0) annotation(
        Placement(visible = true, transformation(origin = {-6, 34}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant const_v(k = 0) annotation(
        Placement(visible = true, transformation(origin = {-94, 84}, extent = {{-7, 7}, {7, -7}}, rotation = -180)));
      Modelica.Blocks.Sources.Constant const_d(k = 0) annotation(
        Placement(visible = true, transformation(origin = {-94, 64}, extent = {{-7, 7}, {7, -7}}, rotation = -180)));
      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a HeatPortCon annotation(
        Placement(visible = true, transformation(origin = {-182, 4}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-2, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a HeatPortRad annotation(
        Placement(visible = true, transformation(origin = {-182, -22}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-2, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.HeatTransfer.Celsius.TemperatureSensor temperatureSensor annotation(
        Placement(visible = true, transformation(origin = {58, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealOutput T annotation(
        Placement(visible = true, transformation(origin = {98, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {150, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      // Lossnay
      Buildings.Fluid.HeatExchangers.ConstantEffectiveness hex(redeclare package Medium1 = Medium, redeclare package Medium2 = Medium, dp1_nominal = 100, dp2_nominal = 100, eps = 0.7, m1_flow_nominal = 2 * (room_w * room_d * room_h) * 1.2 * 0.3 / 3600, m2_flow_nominal = 2 * (room_w * room_d * room_h) * 1.2 * 0.3 / 3600) "Lossnay" annotation(
        Placement(visible = true, transformation(extent = {{-78, -98}, {-58, -78}}, rotation = 0)));
      Buildings.Fluid.Movers.FlowControlled_m_flow fanSup(redeclare package Medium = Medium, allowFlowReversal = true, energyDynamics = Modelica.Fluid.Types.Dynamics.SteadyState, m_flow_nominal = 2 * (room_w * room_d * room_h) * 1.2 * 0.37 / 3600, nominalValuesDefineDefaultPressureCurve = true, use_inputFilter = false) "Supply air fan" annotation(
        Placement(visible = true, transformation(extent = {{-142, -92}, {-122, -72}}, rotation = 0)));
      Buildings.Fluid.Movers.FlowControlled_m_flow fanRet(redeclare package Medium = Medium, T_start = mediumRoomAir_initT, addPowerToMedium = true, allowFlowReversal = true, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, inputType = Buildings.Fluid.Types.InputType.Continuous, m_flow_nominal = m_flow_nominal_lossnay, massFlowRates = {0, 0.5, 1} * m_flow_nominal_lossnay, nominalValuesDefineDefaultPressureCurve = true, use_inputFilter = false) "Return air fan" annotation(
        Placement(visible = true, transformation(origin = {-132, -122}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant const(k = 2 * (room_w * room_d * room_h) * 1.2 * 0.37 / 3600) annotation(
        Placement(visible = true, transformation(origin = {-156, -52}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Fluid.Sources.Boundary_pT out(redeclare package Medium = Medium, T = amb_T, p = amb_P, nPorts = 2) "Outside air conditions" annotation(
        Placement(visible = true, transformation(origin = {-192, -94}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Buildings.Fluid.FixedResistances.PressureDrop dpRet(from_dp = false, redeclare package Medium = Medium, m_flow_nominal = 6 * 4 * 3 * 1.2 * 0.3 / 3600, dp_nominal = 10) "Pressure drop at facade" annotation(
        Placement(visible = true, transformation(origin = {-100, -122}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      Buildings.Fluid.FixedResistances.PressureDrop dpSup(from_dp = false, redeclare package Medium = Medium, m_flow_nominal = 6 * 4 * 3 * 1.2 * 0.3 / 3600, dp_nominal = 10) "Pressure drop at facade" annotation(
        Placement(visible = true, transformation(origin = {-100, -82}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(thermalConductor_wall1.port_b, heatCapacitor.port) annotation(
        Line(points = {{-82, 34}, {-74, 34}, {-74, 44}, {-68, 44}}, color = {191, 0, 0}));
      connect(heatCapacitor.port, thermalConductor_wall2.port_b) annotation(
        Line(points = {{-68, 44}, {-62, 44}, {-62, 34}, {-54, 34}}, color = {191, 0, 0}));
      connect(con_ext.solid, thermalConductor_wall1.port_a) annotation(
        Line(points = {{-122, 34}, {-102, 34}}, color = {191, 0, 0}));
      connect(con_ext.fluid, ambT.port) annotation(
        Line(points = {{-142, 34}, {-152, 34}}, color = {191, 0, 0}));
      connect(thermalConductor_wall2.port_a, con_int.solid) annotation(
        Line(points = {{-34, 34}, {-16, 34}}, color = {191, 0, 0}));
      connect(con_int.fluid, roomAir.heatPort) annotation(
        Line(points = {{4, 34}, {8, 34}, {8, 2}, {12, 2}}, color = {191, 0, 0}));
      connect(const_v.y, con_ext.v) annotation(
        Line(points = {{-101.7, 84}, {-111.7, 84}, {-111.7, 44}, {-119.7, 44}}, color = {0, 0, 127}));
      connect(const_d.y, con_ext.dir) annotation(
        Line(points = {{-101.7, 64}, {-109.7, 64}, {-109.7, 40}, {-119.7, 40}}, color = {0, 0, 127}));
      connect(HeatPortCon, roomAir.heatPort) annotation(
        Line(points = {{-182, 4}, {12, 4}, {12, 2}}, color = {191, 0, 0}));
      connect(HeatPortRad, roomAir.heatPort) annotation(
        Line(points = {{-182, -22}, {-10, -22}, {-10, 2}, {12, 2}}, color = {191, 0, 0}));
      connect(roomAir.heatPort, temperatureSensor.port) annotation(
        Line(points = {{12, 2}, {4, 2}, {4, -20}, {48, -20}, {48, 2}}, color = {191, 0, 0}));
      connect(T, temperatureSensor.T) annotation(
        Line(points = {{98, 2}, {68, 2}}, color = {0, 0, 127}));
      connect(roomAir.ports[1], hex.port_a2) annotation(
        Line(points = {{22, -8}, {22, -94}, {-58, -94}}, color = {0, 127, 255}));
      connect(roomAir.ports[2], hex.port_b1) annotation(
        Line(points = {{22, -8}, {20, -8}, {20, -82}, {-58, -82}}, color = {0, 127, 255}));
      connect(const.y, fanSup.m_flow_in) annotation(
        Line(points = {{-145, -52}, {-132, -52}, {-132, -70}}, color = {0, 0, 127}));
      connect(const.y, fanRet.m_flow_in) annotation(
        Line(points = {{-145, -52}, {-132, -52}, {-132, -110}}, color = {0, 0, 127}));
      connect(out.ports[1], fanSup.port_a) annotation(
        Line(points = {{-182, -94}, {-166, -94}, {-166, -82}, {-142, -82}}, color = {0, 127, 255}));
      connect(fanSup.port_b, dpSup.port_a) annotation(
        Line(points = {{-122, -82}, {-110, -82}}, color = {0, 127, 255}));
      connect(dpSup.port_b, hex.port_a1) annotation(
        Line(points = {{-90, -82}, {-78, -82}}, color = {0, 127, 255}));
      connect(out.ports[2], fanRet.port_b) annotation(
        Line(points = {{-182, -94}, {-166, -94}, {-166, -122}, {-142, -122}}, color = {0, 127, 255}));
      connect(fanRet.port_a, dpRet.port_b) annotation(
        Line(points = {{-122, -122}, {-110, -122}}, color = {0, 127, 255}));
      connect(dpRet.port_a, hex.port_b2) annotation(
        Line(points = {{-90, -122}, {-84, -122}, {-84, -94}, {-78, -94}}, color = {0, 127, 255}));
      annotation(
        Icon(graphics = {Rectangle(lineColor = {0, 85, 255}, fillColor = {0, 170, 255}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-100, 100}, {100, -100}}), Text(origin = {2, 133}, lineColor = {0, 0, 255}, extent = {{-104, 21}, {104, -21}}, textString = "%name"), Line(origin = {101.574, -38}, points = {{-15, 0}, {15, 0}}, color = {0, 0, 255}, thickness = 0.75, arrow = {Arrow.None, Arrow.Filled}, arrowSize = 7), Line(origin = {71.2439, -63.1815}, points = {{45, 0}, {15, 0}}, color = {0, 0, 255}, thickness = 0.75, arrow = {Arrow.None, Arrow.Filled}, arrowSize = 7), Rectangle(origin = {0, 110}, fillColor = {103, 103, 103}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, lineThickness = 0, extent = {{-100, 10}, {100, -10}}), Rectangle(origin = {0, -110}, fillColor = {103, 103, 103}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, lineThickness = 0, extent = {{-100, 10}, {100, -10}}), Rectangle(origin = {-110, 0}, fillColor = {103, 103, 103}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, lineThickness = 0, extent = {{-10, 120}, {10, -120}}), Rectangle(origin = {110, 50}, fillColor = {103, 103, 103}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, lineThickness = 0, extent = {{-10, 70}, {10, -70}}), Line(origin = {117.082, -0.725368}, points = {{-23, 0}, {23, 0}, {23, 0}}), Rectangle(origin = {110, -100}, fillColor = {103, 103, 103}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, lineThickness = 0, extent = {{-10, 20}, {10, -20}})}, coordinateSystem(extent = {{-120, 160}, {160, -120}})),
        Diagram(coordinateSystem(extent = {{-220, 100}, {120, -100}})));
    end RoomWithLossnay;

    model RoomCycle1ZoneTempCtrl
      replaceable package Medium = Buildings.Media.Water;
      inner Modelica.Fluid.System system annotation(
        Placement(visible = true, transformation(origin = {-172, 110}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      // Parameters
      parameter Modelica.Units.SI.Temperature amb_T = 273.15 + 5 "Ambience initial temperature";
      // Initialization
      parameter Modelica.Units.SI.Temperature medium_initT = 273.15 + 30 "Water initial temperature";
      parameter Modelica.Units.SI.Temperature rad_initT = 273.15 + 40 "Radiator water initial temperature";
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
      //
      final parameter Modelica.Units.SI.MassFlowRate m_flow_nominal = 0.1 "Mass flow rate";
      final parameter Modelica.Units.SI.PressureDifference dp_nominal = 4500 "Design pressure drop";
      // Components
      Buildings.Fluid.Movers.FlowControlled_m_flow pump_m_flow_r1(redeclare package Medium = Medium, T_start = medium_initT, allowFlowReversal = false, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, inputType = Buildings.Fluid.Types.InputType.Continuous, m_flow_nominal = m_flow_nominal, massFlowRates = {0, 0.5, 1} * m_flow_nominal, nominalValuesDefineDefaultPressureCurve = true, use_inputFilter = false) "Pump with m_flow input" annotation(
        Placement(visible = true, transformation(origin = {-126, 42}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Buildings.Fluid.HeatExchangers.Radiators.RadiatorEN442_2 rad(redeclare package Medium = Medium, Q_flow_nominal = q_flow_nominal, TAir_nominal = 293.15, T_a_nominal = 353.15, T_b_nominal = 333.15, T_start = rad_initT, allowFlowReversal = false, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, fraRad = 0, p_start = 101325, nEle = 1) annotation(
        Placement(visible = true, transformation(origin = {-26, 10}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      Buildings.Fluid.FixedResistances.PlugFlowPipe pipe_r1(redeclare package Medium = Medium, T_start_in = medium_initT, T_start_out = medium_initT, allowFlowReversal = false, cPip = pip_c, dIns = pip_thickness, dh = pip_dh, kIns = pip_kIns, length = pip_len, rhoPip = pip_rho, thickness = pip_thickness, m_flow_nominal = m_flow_nominal, R = Modelica.Math.log((pip_dh + pip_thickness) / pip_dh) / (2 * Modelica.Constants.pi * pip_kIns) + 1 / pip_air_htc / Modelica.Constants.pi / (pip_dh + pip_thickness)) annotation(
        Placement(visible = true, transformation(origin = {-154, 42}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Buildings.Fluid.FixedResistances.PlugFlowPipe pipe_r2(redeclare package Medium = Medium, T_start_in = medium_initT, T_start_out = medium_initT, allowFlowReversal = false, cPip = pip_c, dIns = pip_thickness, dh = pip_dh, kIns = pip_kIns, length = pip_len, rhoPip = pip_rho, thickness = pip_thickness, m_flow_nominal = m_flow_nominal, R = Modelica.Math.log((pip_dh + pip_thickness) / pip_dh) / (2 * Modelica.Constants.pi * pip_kIns) + 1 / pip_air_htc / Modelica.Constants.pi / (pip_dh + pip_thickness)) annotation(
        Placement(visible = true, transformation(origin = {-52, 42}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Buildings.Fluid.FixedResistances.PlugFlowPipe pipe_r3(redeclare package Medium = Medium, T_start_in = medium_initT, T_start_out = medium_initT, allowFlowReversal = false, cPip = pip_c, dIns = pip_thickness, dh = pip_dh, kIns = pip_kIns, length = pip_len, rhoPip = pip_rho, thickness = pip_thickness, m_flow_nominal = m_flow_nominal, R = Modelica.Math.log((pip_dh + pip_thickness) / pip_dh) / (2 * Modelica.Constants.pi * pip_kIns) + 1 / pip_air_htc / Modelica.Constants.pi / (pip_dh + pip_thickness)) annotation(
        Placement(visible = true, transformation(origin = {-52, -40}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      Buildings.Fluid.FixedResistances.PlugFlowPipe pipe_r4(redeclare package Medium = Medium, T_start_in = medium_initT, T_start_out = medium_initT, allowFlowReversal = false, cPip = pip_c, dIns = pip_thickness, dh = pip_dh, kIns = pip_kIns, length = pip_len, rhoPip = pip_rho, thickness = pip_thickness, m_flow_nominal = m_flow_nominal, R = Modelica.Math.log((pip_dh + pip_thickness) / pip_dh) / (2 * Modelica.Constants.pi * pip_kIns) + 1 / pip_air_htc / Modelica.Constants.pi / (pip_dh + pip_thickness)) annotation(
        Placement(visible = true, transformation(origin = {-154, -40}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      Buildings.Fluid.FixedResistances.PlugFlowPipe pipe_r5(redeclare package Medium = Medium, R = Modelica.Math.log((pip_dh + pip_thickness) / pip_dh) / (2 * Modelica.Constants.pi * pip_kIns) + 1 / pip_air_htc / Modelica.Constants.pi / (pip_dh + pip_thickness), T_start_in = medium_initT, T_start_out = medium_initT, allowFlowReversal = true, cPip = pip_c, dIns = pip_thickness, dh = pip_dh, kIns = pip_kIns, length = pip5_len, m_flow_nominal = m_flow_nominal, rhoPip = pip_rho, thickness = pip_thickness) annotation(
        Placement(visible = true, transformation(origin = {-94, -14}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
      Buildings.Fluid.FixedResistances.Junction jun_r1(redeclare package Medium = Medium, dp_nominal = {0, 0, 0}, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, m_flow_nominal = m_flow_nominal * {1, -1, 1}) annotation(
        Placement(visible = true, transformation(origin = {-94, 42}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Buildings.Fluid.FixedResistances.Junction jun_r2(redeclare package Medium = Medium, dp_nominal = {0, 0, 0}, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, m_flow_nominal = m_flow_nominal * {1, -1, -1}) annotation(
        Placement(visible = true, transformation(origin = {-94, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
      Buildings.HeatTransfer.Sources.FixedTemperature outdoor(T = amb_T) "Boundary temperature" annotation(
        Placement(visible = true, transformation(origin = {-190, 85}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant const_r1(k = 0) annotation(
        Placement(visible = true, transformation(origin = {-68, 3}, extent = {{-7, 7}, {7, -7}}, rotation = 90)));
      Buildings.Fluid.Actuators.Valves.TwoWayLinear val_r1(redeclare package Medium = Medium, CvData = Buildings.Fluid.Types.CvTypes.OpPoint, dpValve_nominal(displayUnit = "kPa") = dp_nominal, m_flow_nominal = m_flow_nominal, use_inputFilter = false) annotation(
        Placement(visible = true, transformation(origin = {-94, 14}, extent = {{10, -10}, {-10, 10}}, rotation = -90)));
      // Interfaces
      //
      Modelica.Fluid.Interfaces.FluidPort_a port_a(redeclare package Medium = Medium) annotation(
        Placement(visible = true, transformation(origin = {-180, 42}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Fluid.Interfaces.FluidPort_b port_b(redeclare package Medium = Medium) annotation(
        Placement(visible = true, transformation(origin = {-180, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealInput m_flow_rate annotation(
        Placement(visible = true, transformation(origin = {-222, 70}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-120, 100}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a port_rad annotation(
        Placement(visible = true, transformation(origin = {24, -2}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a port_con annotation(
        Placement(visible = true, transformation(origin = {24, 22}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(pipe_r1.port_b, pump_m_flow_r1.port_a) annotation(
        Line(points = {{-144, 42}, {-136, 42}}, color = {0, 127, 255}));
      connect(pump_m_flow_r1.port_b, jun_r1.port_1) annotation(
        Line(points = {{-116, 42}, {-104, 42}}, color = {0, 127, 255}));
      connect(jun_r1.port_2, pipe_r2.port_a) annotation(
        Line(points = {{-84, 42}, {-62, 42}}, color = {0, 127, 255}));
      connect(pipe_r2.port_b, rad.port_a) annotation(
        Line(points = {{-42, 42}, {-26, 42}, {-26, 20}}, color = {0, 127, 255}));
      connect(rad.port_b, pipe_r3.port_a) annotation(
        Line(points = {{-26, 0}, {-26, -40}, {-42, -40}}, color = {0, 127, 255}));
      connect(jun_r1.port_3, val_r1.port_b) annotation(
        Line(points = {{-94, 32}, {-94, 24}}, color = {0, 127, 255}));
      connect(val_r1.port_a, pipe_r5.port_b) annotation(
        Line(points = {{-94, 4}, {-94, -4}}, color = {0, 127, 255}));
      connect(pipe_r5.port_a, jun_r2.port_3) annotation(
        Line(points = {{-94, -24}, {-94, -30}}, color = {0, 127, 255}));
      connect(jun_r2.port_2, pipe_r4.port_a) annotation(
        Line(points = {{-104, -40}, {-144, -40}}, color = {0, 127, 255}));
      connect(pipe_r3.port_b, jun_r2.port_1) annotation(
        Line(points = {{-62, -40}, {-84, -40}}, color = {0, 127, 255}));
      connect(pipe_r4.port_b, port_b) annotation(
        Line(points = {{-164, -40}, {-180, -40}}, color = {0, 127, 255}));
      connect(port_a, pipe_r1.port_a) annotation(
        Line(points = {{-180, 42}, {-164, 42}}));
      connect(m_flow_rate, pump_m_flow_r1.m_flow_in) annotation(
        Line(points = {{-222, 70}, {-126, 70}, {-126, 54}}, color = {0, 0, 127}));
      connect(const_r1.y, val_r1.y) annotation(
        Line(points = {{-68, 10}, {-82, 10}, {-82, 14}}, color = {0, 0, 127}));
      connect(rad.heatPortCon, port_con) annotation(
        Line(points = {{-18, 12}, {12, 12}, {12, 22}, {24, 22}}, color = {191, 0, 0}));
      connect(rad.heatPortRad, port_rad) annotation(
        Line(points = {{-18, 8}, {12, 8}, {12, -2}, {24, -2}}, color = {191, 0, 0}));
      connect(outdoor.port, pipe_r1.heatPort) annotation(
        Line(points = {{-182, 86}, {-154, 86}, {-154, 52}}, color = {191, 0, 0}));
      connect(outdoor.port, pipe_r2.heatPort) annotation(
        Line(points = {{-182, 86}, {-52, 86}, {-52, 52}}, color = {191, 0, 0}));
      connect(outdoor.port, pipe_r4.heatPort) annotation(
        Line(points = {{-182, 86}, {-154, 86}, {-154, -30}}, color = {191, 0, 0}));
      connect(outdoor.port, pipe_r3.heatPort) annotation(
        Line(points = {{-182, 86}, {-52, 86}, {-52, -30}}, color = {191, 0, 0}));
      connect(pipe_r5.heatPort, outdoor.port) annotation(
        Line(points = {{-104, -14}, {-154, -14}, {-154, 86}, {-182, 86}}, color = {191, 0, 0}));
      annotation(
        Icon(graphics = {Rectangle(origin = {0, -2}, lineColor = {90, 90, 90}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-100, 100}, {100, -100}}), Text(origin = {0, 120}, lineColor = {0, 0, 255}, extent = {{-100, 20}, {100, -20}}, textString = "%name"), Line(origin = {-25, 25}, points = {{-65, 25}, {65, 25}, {65, -25}, {65, -25}, {65, -25}, {65, -25}}, color = {255, 0, 0}, thickness = 0.5), Line(origin = {-15, 35}, points = {{-75, -85}, {55, -85}, {55, -35}}, color = {0, 0, 255}, thickness = 0.5), Rectangle(extent = {{-38, 100}, {-38, 100}}), Rectangle(origin = {40, 0}, fillColor = {99, 99, 99}, fillPattern = FillPattern.Solid, extent = {{-10, 22}, {10, -22}}), Line(origin = {-32, 27}, points = {{0, 23}, {0, -23}, {0, -23}}, color = {255, 0, 0}, thickness = 0.5), Line(origin = {-31.62, -25.26}, points = {{0, 23}, {0, -23}, {0, -23}}, color = {0, 0, 255}, thickness = 0.5), Polygon(rotation = 90, fillColor = DynamicSelect({0, 0, 0}, {y * 255, y * 255, y * 255}), fillPattern = FillPattern.Solid, points = {{6, 32}, {-2, 38}, {-2, 26}, {6, 32}}), Polygon(rotation = 90, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, points = {{6, 32}, {14, 38}, {14, 26}, {6, 32}})}, coordinateSystem(extent = {{-100, -100}, {100, 100}})),
        Diagram(coordinateSystem(extent = {{-240, 120}, {40, -60}})));
    end RoomCycle1ZoneTempCtrl;

    model RoomCycle1ZoneTempCtrl2
      replaceable package Medium = Buildings.Media.Water;
      inner Modelica.Fluid.System system annotation(
        Placement(visible = true, transformation(origin = {-172, 110}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      // Parameters
      parameter Modelica.Units.SI.Temperature amb_T = 273.15 + 5 "Ambience initial temperature";
      // Initialization
      parameter Modelica.Units.SI.Temperature medium_initT = 273.15 + 30 "Water initial temperature";
      parameter Modelica.Units.SI.Temperature rad_initT = 273.15 + 40 "Radiator water initial temperature";
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
      //
      final parameter Modelica.Units.SI.MassFlowRate m_flow_nominal = 0.1 "Mass flow rate";
      final parameter Modelica.Units.SI.PressureDifference dp_nominal = 4500 "Design pressure drop";
      // Components
      Buildings.Fluid.Movers.FlowControlled_m_flow pump_m_flow_r1(redeclare package Medium = Medium, T_start = medium_initT, allowFlowReversal = false, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, inputType = Buildings.Fluid.Types.InputType.Continuous, m_flow_nominal = m_flow_nominal, massFlowRates = {0, 0.5, 1} * m_flow_nominal, nominalValuesDefineDefaultPressureCurve = true, use_inputFilter = false) "Pump with m_flow input" annotation(
        Placement(visible = true, transformation(origin = {-126, 42}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Buildings.Fluid.HeatExchangers.Radiators.RadiatorEN442_2 rad(redeclare package Medium = Medium, Q_flow_nominal = q_flow_nominal, TAir_nominal = 293.15, T_a_nominal = 353.15, T_b_nominal = 333.15, T_start = rad_initT, allowFlowReversal = false, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, fraRad = 0, p_start = 101325, nEle = 1) annotation(
        Placement(visible = true, transformation(origin = {-26, 10}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
      Buildings.Fluid.FixedResistances.PlugFlowPipe pipe_r1(redeclare package Medium = Medium, T_start_in = medium_initT, T_start_out = medium_initT, allowFlowReversal = false, cPip = pip_c, dIns = pip_thickness, dh = pip_dh, kIns = pip_kIns, length = pip_len, rhoPip = pip_rho, thickness = pip_thickness, m_flow_nominal = m_flow_nominal, R = Modelica.Math.log((pip_dh + pip_thickness) / pip_dh) / (2 * Modelica.Constants.pi * pip_kIns) + 1 / pip_air_htc / Modelica.Constants.pi / (pip_dh + pip_thickness)) annotation(
        Placement(visible = true, transformation(origin = {-154, 42}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Buildings.Fluid.FixedResistances.PlugFlowPipe pipe_r2(redeclare package Medium = Medium, T_start_in = medium_initT, T_start_out = medium_initT, allowFlowReversal = false, cPip = pip_c, dIns = pip_thickness, dh = pip_dh, kIns = pip_kIns, length = pip_len, rhoPip = pip_rho, thickness = pip_thickness, m_flow_nominal = m_flow_nominal, R = Modelica.Math.log((pip_dh + pip_thickness) / pip_dh) / (2 * Modelica.Constants.pi * pip_kIns) + 1 / pip_air_htc / Modelica.Constants.pi / (pip_dh + pip_thickness)) annotation(
        Placement(visible = true, transformation(origin = {-52, 42}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Buildings.Fluid.FixedResistances.PlugFlowPipe pipe_r3(redeclare package Medium = Medium, T_start_in = medium_initT, T_start_out = medium_initT, allowFlowReversal = false, cPip = pip_c, dIns = pip_thickness, dh = pip_dh, kIns = pip_kIns, length = pip_len, rhoPip = pip_rho, thickness = pip_thickness, m_flow_nominal = m_flow_nominal, R = Modelica.Math.log((pip_dh + pip_thickness) / pip_dh) / (2 * Modelica.Constants.pi * pip_kIns) + 1 / pip_air_htc / Modelica.Constants.pi / (pip_dh + pip_thickness)) annotation(
        Placement(visible = true, transformation(origin = {-52, -40}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      Buildings.Fluid.FixedResistances.PlugFlowPipe pipe_r4(redeclare package Medium = Medium, T_start_in = medium_initT, T_start_out = medium_initT, allowFlowReversal = false, cPip = pip_c, dIns = pip_thickness, dh = pip_dh, kIns = pip_kIns, length = pip_len, rhoPip = pip_rho, thickness = pip_thickness, m_flow_nominal = m_flow_nominal, R = Modelica.Math.log((pip_dh + pip_thickness) / pip_dh) / (2 * Modelica.Constants.pi * pip_kIns) + 1 / pip_air_htc / Modelica.Constants.pi / (pip_dh + pip_thickness)) annotation(
        Placement(visible = true, transformation(origin = {-154, -40}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      Buildings.HeatTransfer.Sources.FixedTemperature outdoor(T = amb_T) "Boundary temperature" annotation(
        Placement(visible = true, transformation(origin = {-190, 85}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
      // Interfaces
      //
      Modelica.Fluid.Interfaces.FluidPort_a port_a(redeclare package Medium = Medium) annotation(
        Placement(visible = true, transformation(origin = {-180, 42}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Fluid.Interfaces.FluidPort_b port_b(redeclare package Medium = Medium) annotation(
        Placement(visible = true, transformation(origin = {-180, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Interfaces.RealInput m_flow_rate annotation(
        Placement(visible = true, transformation(origin = {-222, 70}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-120, 100}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a port_rad annotation(
        Placement(visible = true, transformation(origin = {24, -2}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a port_con annotation(
        Placement(visible = true, transformation(origin = {24, 22}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(pipe_r1.port_b, pump_m_flow_r1.port_a) annotation(
        Line(points = {{-144, 42}, {-136, 42}}, color = {0, 127, 255}));
      connect(pipe_r2.port_b, rad.port_a) annotation(
        Line(points = {{-42, 42}, {-26, 42}, {-26, 20}}, color = {0, 127, 255}));
      connect(rad.port_b, pipe_r3.port_a) annotation(
        Line(points = {{-26, 0}, {-26, -40}, {-42, -40}}, color = {0, 127, 255}));
      connect(pipe_r4.port_b, port_b) annotation(
        Line(points = {{-164, -40}, {-180, -40}}, color = {0, 127, 255}));
      connect(port_a, pipe_r1.port_a) annotation(
        Line(points = {{-180, 42}, {-164, 42}}));
      connect(m_flow_rate, pump_m_flow_r1.m_flow_in) annotation(
        Line(points = {{-222, 70}, {-126, 70}, {-126, 54}}, color = {0, 0, 127}));
      connect(rad.heatPortCon, port_con) annotation(
        Line(points = {{-18, 12}, {12, 12}, {12, 22}, {24, 22}}, color = {191, 0, 0}));
      connect(rad.heatPortRad, port_rad) annotation(
        Line(points = {{-18, 8}, {12, 8}, {12, -2}, {24, -2}}, color = {191, 0, 0}));
      connect(outdoor.port, pipe_r1.heatPort) annotation(
        Line(points = {{-182, 86}, {-154, 86}, {-154, 52}}, color = {191, 0, 0}));
      connect(outdoor.port, pipe_r2.heatPort) annotation(
        Line(points = {{-182, 86}, {-52, 86}, {-52, 52}}, color = {191, 0, 0}));
      connect(outdoor.port, pipe_r4.heatPort) annotation(
        Line(points = {{-182, 86}, {-154, 86}, {-154, -30}}, color = {191, 0, 0}));
      connect(outdoor.port, pipe_r3.heatPort) annotation(
        Line(points = {{-182, 86}, {-52, 86}, {-52, -30}}, color = {191, 0, 0}));
      connect(pump_m_flow_r1.port_b, pipe_r2.port_a) annotation(
        Line(points = {{-116, 42}, {-62, 42}}, color = {0, 127, 255}));
      connect(pipe_r3.port_b, pipe_r4.port_a) annotation(
        Line(points = {{-62, -40}, {-144, -40}}, color = {0, 127, 255}));
      annotation(
        Icon(graphics = {Rectangle(lineColor = {90, 90, 90}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, lineThickness = 1, extent = {{-100, 100}, {100, -100}}), Text(origin = {0, 120}, lineColor = {0, 0, 255}, extent = {{-100, 20}, {100, -20}}, textString = "%name"), Line(origin = {-25, 25}, points = {{-65, 25}, {65, 25}, {65, -25}, {65, -25}, {65, -25}, {65, -25}}, color = {255, 0, 0}, thickness = 0.5), Line(origin = {-15, 35}, points = {{-75, -85}, {55, -85}, {55, -35}}, color = {0, 0, 255}, thickness = 0.5), Rectangle(extent = {{-38, 100}, {-38, 100}}), Rectangle(origin = {40, 0}, fillColor = {99, 99, 99}, fillPattern = FillPattern.Solid, extent = {{-10, 22}, {10, -22}})}, coordinateSystem(extent = {{-100, -100}, {100, 100}})),
        Diagram(coordinateSystem(extent = {{-240, 120}, {40, -60}})),
        Documentation(info = "<html><p>
    Air To Waterのうち，部屋側の水回路のモデル．</br>
    1台のポンプと1台のラジエーターを含む．
    </p></html>", revisions = "<html>
    <ul>
    <li>
    </li>
    <li>
    </li>
    <li>
    </li>
    </ul>
    </html>"));
    end RoomCycle1ZoneTempCtrl2;

    model SimpleRoom
      replaceable package Medium = Buildings.Media.Air;
      inner Modelica.Fluid.System system annotation(
        Placement(visible = true, transformation(origin = {-170, 92}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      extends BaseClasses.PartialSimpleAirHeatBalance;
      // Parameters
      parameter Modelica.Units.SI.Temperature amb_T = 273.15 + 5 "Ambience initial temperature";
      // Initialization
      parameter Modelica.Units.SI.Temperature mediumRoomAir_initT = 273.15 + 20 "Room air initial temperature";
      parameter Modelica.Units.SI.Pressure mediumRoomAir_initP = 101325 "Room air initial Pressure";
      parameter Modelica.Units.SI.Temperature rw_initT = 273.15 + 15 "Room wall initial temperature";
      // Room
      parameter Modelica.Units.SI.Length room_w = 1.6 "Width of room";
      parameter Modelica.Units.SI.Length room_d = 4.8 "Depth of room";
      parameter Modelica.Units.SI.Length room_h = 2.4 "Height of room";
      parameter Modelica.Units.SI.Area room_area = 2 * (room_w * room_d + room_d * room_h + room_h * room_w) "Area of room";
      parameter Modelica.Units.SI.Length rw_thickness = 0.1 "Thickness of room wall";
      parameter Modelica.Units.SI.Density rw_rho = 144 "Density of room wall";
      parameter Modelica.Units.SI.SpecificHeatCapacity rw_c = 1168 "Specific heat capacity of room wall";
      parameter Modelica.Units.SI.ThermalConductivity rw_k = 0.3 "Heat conductivity of room wall";
      parameter Modelica.Units.SI.Volume rw_vol = room_area * rw_thickness;
      parameter Modelica.Units.SI.CoefficientOfHeatTransfer rw_air_htc = 5 "Air-Room wall heat transfer coeff.";
      final parameter Modelica.Units.SI.MassFlowRate m_flow_nominal = 0.1 "Mass flow rate";
      //final parameter Modelica.Units.SI.PressureDifference dp_nominal = 4500 "Design pressure drop";
      // Components
      //
      Buildings.HeatTransfer.Sources.FixedTemperature ambT(T = amb_T) annotation(
        Placement(visible = true, transformation(origin = {-176, 46}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
      Buildings.Fluid.MixingVolumes.MixingVolume roomAir(redeclare package Medium = Medium, T_start = mediumRoomAir_initT, V = room_w * room_d * room_h, allowFlowReversal = true, m_flow_nominal = 0.01, massDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, nPorts = 0, p_start(displayUnit = "Pa") = mediumRoomAir_initP, use_C_flow = false) annotation(
        Placement(visible = true, transformation(origin = {22, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.HeatTransfer.Components.HeatCapacitor heatCapacitor(C = rw_vol * rw_rho * rw_c, T(fixed = true, start = rw_initT)) annotation(
        Placement(visible = true, transformation(origin = {-84, 66}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor_wall1(G = rw_k * room_area * 2 / rw_thickness) annotation(
        Placement(visible = true, transformation(origin = {-108, 46}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor_wall2(G = rw_k * room_area * 2 / rw_thickness) annotation(
        Placement(visible = true, transformation(origin = {-60, 46}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      Buildings.HeatTransfer.Convection.Exterior con_ext(A = room_area, azi = 3.14 / 2, conMod = Buildings.HeatTransfer.Types.ExteriorConvection.Fixed, hFixed = rw_air_htc, til = 0) annotation(
        Placement(visible = true, transformation(origin = {-148, 46}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
      Buildings.HeatTransfer.Convection.Interior con_int(A = room_area, conMod = Buildings.HeatTransfer.Types.InteriorConvection.Fixed, hFixed = rw_air_htc, til = 0) annotation(
        Placement(visible = true, transformation(origin = {-22, 46}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Constant const_v(k = 0) annotation(
        Placement(visible = true, transformation(origin = {-110, 96}, extent = {{-7, 7}, {7, -7}}, rotation = -180)));
      Modelica.Blocks.Sources.Constant const_d(k = 0) annotation(
        Placement(visible = true, transformation(origin = {-110, 76}, extent = {{-7, 7}, {7, -7}}, rotation = -180)));
      Modelica.Thermal.HeatTransfer.Celsius.TemperatureSensor temperatureSensor annotation(
        Placement(visible = true, transformation(origin = {48, -22}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
  connect(thermalConductor_wall1.port_b, heatCapacitor.port) annotation(
        Line(points = {{-98, 46}, {-90, 46}, {-90, 56}, {-84, 56}}, color = {191, 0, 0}));
  connect(heatCapacitor.port, thermalConductor_wall2.port_b) annotation(
        Line(points = {{-84, 56}, {-78, 56}, {-78, 46}, {-70, 46}}, color = {191, 0, 0}));
  connect(con_ext.solid, thermalConductor_wall1.port_a) annotation(
        Line(points = {{-138, 46}, {-118, 46}}, color = {191, 0, 0}));
  connect(con_ext.fluid, ambT.port) annotation(
        Line(points = {{-158, 46}, {-168, 46}}, color = {191, 0, 0}));
  connect(thermalConductor_wall2.port_a, con_int.solid) annotation(
        Line(points = {{-50, 46}, {-32, 46}}, color = {191, 0, 0}));
  connect(con_int.fluid, roomAir.heatPort) annotation(
        Line(points = {{-12, 46}, {-12, 2}, {12, 2}}, color = {191, 0, 0}));
  connect(const_v.y, con_ext.v) annotation(
        Line(points = {{-117.7, 96}, {-127.7, 96}, {-127.7, 56}, {-135.7, 56}}, color = {0, 0, 127}));
  connect(const_d.y, con_ext.dir) annotation(
        Line(points = {{-117.7, 76}, {-125.7, 76}, {-125.7, 52}, {-135.7, 52}}, color = {0, 0, 127}));
  connect(roomAir.heatPort, temperatureSensor.port) annotation(
        Line(points = {{12, 2}, {4, 2}, {4, -22}, {38, -22}}, color = {191, 0, 0}));
/*
      for i in 1:nPorts loop
        connect(ports[i],roomAir.ports[i])
        annotation (Line(
          points={{-260,-60},{-218,-60},{-218,-206},{52,-206},{52,-141.9}},
          color={0,127,255},
          smooth=Smooth.None));
      end for;
      */

      connect(HeatPortCon, roomAir.heatPort) annotation(
        Line(points = {{-100, 20}, {-12, 20}, {-12, 2}, {12, 2}}, color = {191, 0, 0}));
  connect(HeatPortRad, roomAir.heatPort) annotation(
        Line(points = {{-100, -20}, {-12, -20}, {-12, 2}, {12, 2}}, color = {191, 0, 0}));
  connect(temperatureSensor.T, T) annotation(
        Line(points = {{58, -22}, {58, 0}, {102, 0}}, color = {0, 0, 127}));
      protected
  annotation(
        Icon(graphics = {Rectangle(origin = {0, 90}, fillColor = {103, 103, 103}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, lineThickness = 0, extent = {{-80, 10}, {80, -10}}), Rectangle(origin = {0, -90}, fillColor = {103, 103, 103}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, lineThickness = 0, extent = {{-80, 10}, {80, -10}}), Rectangle(origin = {-90, 0}, fillColor = {103, 103, 103}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, lineThickness = 0, extent = {{-10, 100}, {10, -100}}), Rectangle(origin = {90, 0}, fillColor = {103, 103, 103}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, lineThickness = 0, extent = {{-10, 100}, {10, -100}})}, coordinateSystem(extent = {{-140, 140}, {140, -100}})),
        Diagram(coordinateSystem(extent = {{-200, 120}, {120, -40}})));
    end SimpleRoom;

    package BaseClasses
      partial model PartialSimpleAirHeatBalance
        Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a HeatPortCon annotation(
          Placement(visible = true, transformation(origin = {-100, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-120, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a HeatPortRad annotation(
          Placement(visible = true, transformation(origin = {-100, -20}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-120, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Blocks.Interfaces.RealOutput T annotation(
          Placement(visible = true, transformation(origin = {102, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {130, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        annotation(
          Icon(graphics = {Rectangle(fillColor = {0, 170, 255}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, extent = {{-80, 80}, {80, -80}}), Line(origin = {-73, 50}, points = {{-37, 0}, {0, 0}}, color = {189, 0, 0}), Line(origin = {-73.1, -49.85}, points = {{-37, 0}, {0, 0}}, color = {189, 0, 0}), Line(origin = {92, 0}, points = {{28, 0}, {-28, 0}}), Text(origin = {0, 120}, lineColor = {0, 0, 255}, extent = {{-140, 20}, {140, -20}}, textString = "%name")}));
      end PartialSimpleAirHeatBalance;
    end BaseClasses;
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

    model Pump "Example model of movers using a parameter for setting the stage"
      extends Modelica.Icons.Example;
      package Medium = Buildings.Media.Water;
      parameter Modelica.Units.SI.MassFlowRate m_flow_nominal = 2 "Nominal mass flow rate";
      Modelica.Blocks.Sources.Ramp ramp(duration = 100, height = 0.9, offset = 0.1) "Ramp input for all movers" annotation(
        Placement(transformation(extent = {{-80, 60}, {-60, 80}})));
      Modelica.Blocks.Math.Gain gai_m_flow(k = m_flow_nominal) "Nominal mass flow rate" annotation(
        Placement(transformation(extent = {{-40, 10}, {-20, 30}})));
      Buildings.Fluid.Movers.FlowControlled_m_flow pump_m_flow(redeclare package Medium = Medium, m_flow_nominal = m_flow_nominal, use_inputFilter = false, massFlowRates = {0, 0.5, 1} * m_flow_nominal, inputType = Buildings.Fluid.Types.InputType.Continuous, energyDynamics = Modelica.Fluid.Types.Dynamics.SteadyState, nominalValuesDefineDefaultPressureCurve = true, allowFlowReversal = false) "Pump with m_flow input" annotation(
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

    model HeatPump
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
      System.MyComponents.HeatPump hex_hp(K = 2, redeclare package Medium = Medium, T = 100, InitialTemp = 273.15 + 10) annotation(
        Placement(visible = true, transformation(origin = {16, 6}, extent = {{-12, 14}, {12, 36}}, rotation = 0)));
      Modelica.Blocks.Sources.Ramp ramp1(duration = 0, height = 3000, offset = 0, startTime = 10) annotation(
        Placement(transformation(origin = {-32, 60}, extent = {{-10, -10}, {10, 10}})));
      Buildings.Fluid.FixedResistances.LosslessPipe pip1(redeclare package Medium = Medium, allowFlowReversal = false) annotation(
        Placement(transformation(origin = {52, 0}, extent = {{-10, -10}, {10, 10}})));
    equation
      connect(ramp.y, boundary1.m_flow_in) annotation(
        Line(points = {{-75, 8}, {-56, 8}}, color = {0, 0, 127}));
      connect(boundary1.ports[1], hex_hp.port_a) annotation(
        Line(points = {{-36, 0}, {11, 0}}, color = {0, 127, 255}));
      connect(hex_hp.port_a, temperature_a.port) annotation(
        Line(points = {{11, 0}, {-2, 0}, {-2, 8}, {-22, 8}}, color = {0, 127, 255}));
      connect(ramp1.y, hex_hp.freq) annotation(
        Line(points = {{-21, 60}, {-14, 60}, {-14, 10}, {8, 10}}, color = {0, 0, 127}));
      connect(temperature_b.port, hex_hp.port_b) annotation(
        Line(points = {{38, 10}, {38, 0}, {27, 0}}, color = {0, 127, 255}));
      connect(hex_hp.port_b, pip1.port_a) annotation(
        Line(points = {{27, 0}, {42, 0}}, color = {0, 127, 255}));
      connect(pip1.port_b, boundary.ports[1]) annotation(
        Line(points = {{62, 0}, {72, 0}}, color = {0, 127, 255}));
      annotation(
        experiment(StartTime = 0, StopTime = 150, Tolerance = 1e-6, Interval = 0.3),
        Diagram(coordinateSystem(extent = {{-100, 80}, {80, -20}})));
    end HeatPump;

    model HexControl
      extends Modelica.Icons.Example;
      System.MyComponents.CompressorController hexControl annotation(
        Placement(visible = true, transformation(origin = {-6, 14}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression(y = 30) annotation(
        Placement(visible = true, transformation(origin = {-48, 22}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.Sine sine(amplitude = 35, f = 0.01, startTime = 1) annotation(
        Placement(visible = true, transformation(origin = {-48, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(realExpression.y, hexControl.T_set) annotation(
        Line(points = {{-36, 22}, {-16, 22}, {-16, 20}}, color = {0, 0, 127}));
      connect(sine.y, hexControl.T_BT) annotation(
        Line(points = {{-36, 2}, {-22, 2}, {-22, 8}, {-16, 8}}, color = {0, 0, 127}));
      annotation(
        experiment(StopTime = 100, StartTime = 0, Tolerance = 1e-06, Interval = 0.2));
    end HexControl;

    model RoomCycle1ZoneTempCtrl
      extends Modelica.Icons.Example;
      replaceable package Medium = Buildings.Media.Water;
      System.MyComponents.RoomCycle1ZoneTempCtrl roomCycle1ZoneTempCtrl(redeclare package Medium = Medium) annotation(
        Placement(visible = true, transformation(origin = {37, -21}, extent = {{-19, -19}, {19, 19}}, rotation = 0)));
      Modelica.Fluid.Sources.Boundary_pT boundary(redeclare package Medium = Medium, T = 273.15 + 70, nPorts = 1, p = 101400) annotation(
        Placement(visible = true, transformation(origin = {-44, 6}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Fluid.Sources.Boundary_pT boundary1(redeclare package Medium = Medium, nPorts = 1, p = 101325) annotation(
        Placement(visible = true, transformation(origin = {-36, -28}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Modelica.Blocks.Sources.RealExpression realExpression(y = 0.1) annotation(
        Placement(visible = true, transformation(origin = {-14, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    equation
      connect(boundary1.ports[1], roomCycle1ZoneTempCtrl.port_b) annotation(
        Line(points = {{-26, -28}, {18, -28}, {18, -30}}, color = {0, 127, 255}));
      connect(boundary.ports[1], roomCycle1ZoneTempCtrl.port_a) annotation(
        Line(points = {{-34, 6}, {18, 6}, {18, -12}}, color = {0, 127, 255}));
      connect(realExpression.y, roomCycle1ZoneTempCtrl.m_flow_rate) annotation(
        Line(points = {{-2, 38}, {14, 38}, {14, -2}}, color = {0, 0, 127}));
      annotation(
        experiment(StopTime = 100, StartTime = 0, Tolerance = 1e-06, Interval = 0.2),
        Diagram);
    end RoomCycle1ZoneTempCtrl;

    model RoomWithLossnay
      extends Modelica.Icons.Example;
      replaceable package Medium = Buildings.Media.Air;
      MyComponents.RoomWithLossnay roomWithLossnay(redeclare package Medium = Medium) annotation(
        Placement(visible = true, transformation(origin = {22, -4}, extent = {{-10, 16}, {14, 42}}, rotation = 0)));
    equation

      annotation(
        experiment(StopTime = 100, StartTime = 0, Tolerance = 1e-06, Interval = 0.2),
        Diagram(coordinateSystem(extent = {{0, 40}, {40, -20}})));
    end RoomWithLossnay;
  end Test;

  model plant "Hydronic heating model"
    extends Modelica.Icons.Example;
    replaceable package Medium = Buildings.Media.Water;
    replaceable package MediumAir = Buildings.Media.Air;
    inner Modelica.Fluid.System system annotation(
      Placement(visible = true, transformation(origin = {-142, 140}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    import SI = Modelica.Units.SI;
    // Parameters
    import Modelica.Constants.pi;
    parameter Modelica.Units.SI.Temperature Amb_T = 273.15 + 5 "Ambience temperature";
    // Initialization
    parameter Modelica.Units.SI.Temperature medium_initT = 273.15 + 30 "Water initial temperature";
    parameter Modelica.Units.SI.Temperature BT_initT = 273.15 + 30 "Buffer tank init temperature";
    parameter Modelica.Units.SI.Pressure BT_initP = 101325 "Buffer tank init pressure";
    parameter Modelica.Units.SI.Temperature rad_initT = 273.15 + 40 "Radiator water initial temperature";
    // Buffer Tank
    parameter Modelica.Units.SI.Volume BT_vol = 0.2 "Buffer tank volume";
    parameter Modelica.Units.SI.Length BT_height = 1 "Buffer tank height";
    // Pipe
    parameter Modelica.Units.SI.Length pip_len = 1 "Length of pipe";
    parameter Modelica.Units.SI.Length pipr5_len = 0.2 "Length of pipe wall";
    parameter Modelica.Units.SI.SpecificHeatCapacity pip_c = 386 "Specific heat of pipe wall material";
    parameter Modelica.Units.SI.Density pip_rho = 8960 "Density of pipe wall material";
    parameter Modelica.Units.SI.Length pip_dh = 0.025 "Inner diameter of pipe";
    parameter Modelica.Units.SI.Length pip_thickness = 0.005 "Thickness of pipe wall";
    parameter Modelica.Units.SI.ThermalConductivity pip_kIns = 398 "Heat conductivity of pipe insulation";
    parameter Modelica.Units.SI.CoefficientOfHeatTransfer pip_air_htc = 10 "Air-Pipe heat transfer coefficient";
    // Pump
    parameter Real m_flow = 10 / 60 * 0.001 * 1000 "Mass flow";
    parameter Real m_flow_room = 5 / 60 * 0.001 * 1000 "Mass flow of room cycle";
    // Radiator
    parameter Modelica.Units.SI.Power q_flow_nominal = 5000 "Rated heat dissipation amount";
    //
    final parameter Modelica.Units.SI.MassFlowRate m_flow_nominal = 0.1 "Mass flow rate";
    final parameter Modelica.Units.SI.PressureDifference dp_nominal = 4500 "Design pressure drop";
    //===========
    // components
    // heat pump
    Buildings.Fluid.Movers.FlowControlled_m_flow pump_m_flow(redeclare package Medium = Medium, T_start = medium_initT, allowFlowReversal = false, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, inputType = Buildings.Fluid.Types.InputType.Continuous, m_flow_nominal = m_flow_nominal, massFlowRates = {0, 0.5, 1} * m_flow_nominal, nominalValuesDefineDefaultPressureCurve = true, use_inputFilter = false) "Pump with m_flow input" annotation(
      Placement(transformation(origin = {-54, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
    Buildings.Fluid.Storage.Stratified BufferTank(redeclare package Medium = Medium, T_start = BT_initT, VTan = BT_vol, dIns = 0.01, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, hTan = BT_height, m_flow_nominal = m_flow_nominal, nSeg = 2) "comment" annotation(
      Placement(transformation(origin = {27, 81}, extent = {{15, -15}, {-15, 15}})));
    Buildings.Fluid.FixedResistances.PlugFlowPipe pipe1(redeclare package Medium = Medium, length = pip_len, T_start_in = medium_initT, T_start_out = medium_initT, allowFlowReversal = false, cPip = pip_c, rhoPip = pip_rho, thickness = pip_thickness, dh = pip_dh, dIns = pip_thickness, kIns = pip_kIns, m_flow_nominal = m_flow_nominal, R = Modelica.Math.log((pip_dh + pip_thickness) / pip_dh) / (2 * Modelica.Constants.pi * pip_kIns) + 1 / pip_air_htc / Modelica.Constants.pi / (pip_dh + pip_thickness)) "Outside PlugFlowPipe" annotation(
      Placement(transformation(origin = {-74, 46}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
    Buildings.Fluid.FixedResistances.PlugFlowPipe pipe2(redeclare package Medium = Medium, T_start_in = medium_initT, T_start_out = medium_initT, allowFlowReversal = false, cPip = pip_c, length = pip_len, rhoPip = pip_rho, thickness = pip_thickness, dh = pip_dh, dIns = pip_thickness, kIns = pip_kIns, m_flow_nominal = m_flow_nominal, R = Modelica.Math.log((pip_dh + pip_thickness) / pip_dh) / (2 * Modelica.Constants.pi * pip_kIns) + 1 / pip_air_htc / Modelica.Constants.pi / (pip_dh + pip_thickness)) "Indoor PlugFlowPipe" annotation(
      Placement(transformation(origin = {-24, 96}, extent = {{-10, -10}, {10, 10}})));
    Buildings.Fluid.FixedResistances.PlugFlowPipe pipe3(redeclare package Medium = Medium, T_start_in = medium_initT, T_start_out = medium_initT, allowFlowReversal = false, cPip = pip_c, length = pip_len, rhoPip = pip_rho, thickness = pip_thickness, dh = pip_dh, dIns = pip_thickness, kIns = pip_kIns, m_flow_nominal = m_flow_nominal, R = Modelica.Math.log((pip_dh + pip_thickness) / pip_dh) / (2 * Modelica.Constants.pi * pip_kIns) + 1 / pip_air_htc / Modelica.Constants.pi / (pip_dh + pip_thickness)) "comment" annotation(
      Placement(transformation(origin = {26, 40}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Buildings.Fluid.FixedResistances.PlugFlowPipe pipe4(redeclare package Medium = Medium, T_start_in = medium_initT, T_start_out = medium_initT, allowFlowReversal = false, cPip = pip_c, length = pip_len, rhoPip = pip_rho, thickness = pip_thickness, dh = pip_dh, dIns = pip_thickness, kIns = pip_kIns, m_flow_nominal = m_flow_nominal, R = Modelica.Math.log((pip_dh + pip_thickness) / pip_dh) / (2 * Modelica.Constants.pi * pip_kIns) + 1 / pip_air_htc / Modelica.Constants.pi / (pip_dh + pip_thickness)) "" annotation(
      Placement(transformation(origin = {-74, -14}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
    Buildings.Fluid.FixedResistances.PlugFlowPipe pipe5(redeclare package Medium = Medium, T_start_in = medium_initT, T_start_out = medium_initT, allowFlowReversal = false, cPip = pip_c, length = pip_len, rhoPip = pip_rho, thickness = pip_thickness, dh = pip_dh, dIns = pip_thickness, kIns = pip_kIns, m_flow_nominal = m_flow_nominal, R = Modelica.Math.log((pip_dh + pip_thickness) / pip_dh) / (2 * Modelica.Constants.pi * pip_kIns) + 1 / pip_air_htc / Modelica.Constants.pi / (pip_dh + pip_thickness)) "" annotation(
      Placement(transformation(origin = {-24, 42}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Buildings.Fluid.Actuators.Valves.TwoWayLinear val(redeclare package Medium = Medium, m_flow_nominal = m_flow_nominal, CvData = Buildings.Fluid.Types.CvTypes.OpPoint, dpValve_nominal(displayUnit = "kPa") = dp_nominal, use_inputFilter = false) "Valve model, linear opening characteristics" annotation(
      Placement(visible = true, transformation(origin = {-44, 72}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Constant const(k = 0) annotation(
      Placement(transformation(origin = {-62, 86}, extent = {{-6, -6}, {6, 6}})));
    Buildings.HeatTransfer.Sources.FixedTemperature outdoor(each T = Amb_T) "Boundary temperature" annotation(
      Placement(transformation(origin = {-118, 117}, extent = {{-7, -7}, {7, 7}})));
    Modelica.Blocks.Sources.Ramp ramp1(duration = 0, height = m_flow, offset = 0, startTime = 0) annotation(
      Placement(visible = true, transformation(origin = {-81, -53}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
    Buildings.Fluid.Storage.ExpansionVessel exp(redeclare package Medium = Medium, p_start = 101325, V_start = 0.001, T_start = medium_initT) annotation(
      Placement(transformation(origin = {26, -54}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
    Modelica.Thermal.HeatTransfer.Celsius.TemperatureSensor temperatureSensor1 annotation(
      Placement(visible = true, transformation(origin = {46, 18}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    System.MyComponents.HeatPump heatPump(redeclare package Medium = Medium, InitialTemp = medium_initT) annotation(
      Placement(visible = true, transformation(origin = {-77.6112, 18.3373}, extent = {{-14.6707, 15.2349}, {7.33533, 37.2409}}, rotation = 90)));
    System.MyComponents.CompressorController compCtrl annotation(
      Placement(visible = true, transformation(origin = {-114, -2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.RealExpression Tset(y = 50) annotation(
      Placement(visible = true, transformation(origin = {-146, 4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Buildings.Fluid.FixedResistances.Junction jun(redeclare package Medium = Medium, dp_nominal = {0, 0, 0}, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, m_flow_nominal = m_flow * {1, -1, -1}) "Splitter" annotation(
      Placement(visible = true, transformation(origin = {-74, 72}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
    // romm cycle
    System.MyComponents.RoomWithoutLossnay roomWithoutLossnay1(redeclare package Medium = MediumAir) annotation(
      Placement(visible = true, transformation(origin = {124, 72}, extent = {{-10, 16}, {14, 42}}, rotation = 0)));
    System.MyComponents.PumpController pumpCtrl annotation(
      Placement(visible = true, transformation(origin = {70, 104}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Blocks.Sources.RealExpression TsetAir(y = 25) annotation(
      Placement(visible = true, transformation(origin = {112, 124}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    System.MyComponents.RoomCycle1ZoneTempCtrl2 roomCycle1ZoneTempCtrl1(redeclare package Medium = Medium) annotation(
      Placement(visible = true, transformation(origin = {90, 72}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(pump_m_flow.port_b, pipe4.port_a) annotation(
      Line(points = {{-64, -40}, {-74, -40}, {-74, -24}}, color = {0, 127, 255}));
    connect(val.port_b, pipe5.port_a) annotation(
      Line(points = {{-34, 72}, {-24, 72}, {-24, 52}}, color = {0, 127, 255}));
    connect(pipe5.port_b, pump_m_flow.port_a) annotation(
      Line(points = {{-24, 32}, {-24, -40}, {-44, -40}}, color = {0, 127, 255}));
    connect(const.y, val.y) annotation(
      Line(points = {{-56, 86}, {-44, 86}, {-44, 84}}, color = {0, 0, 127}));
    connect(ramp1.y, pump_m_flow.m_flow_in) annotation(
      Line(points = {{-73, -53}, {-73, -52}, {-54, -52}}, color = {0, 0, 127}));
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
    connect(BufferTank.heaPorVol[1], temperatureSensor1.port) annotation(
      Line(points = {{28, 82}, {47, 82}, {47, 28}, {46, 28}}, color = {191, 0, 0}));
    connect(pipe1.port_a, heatPump.port_b) annotation(
      Line(points = {{-74, 36}, {-74, 32}, {-73, 32}, {-73, 26}}, color = {0, 127, 255}));
    connect(pipe4.port_b, heatPump.port_a) annotation(
      Line(points = {{-74, -4}, {-74, 2}, {-73, 2}, {-73, 11}}, color = {0, 127, 255}));
    connect(compCtrl.y, heatPump.freq) annotation(
      Line(points = {{-104, -2}, {-82, -2}, {-82, 8}}, color = {0, 0, 127}));
    connect(temperatureSensor1.T, compCtrl.T_BT) annotation(
      Line(points = {{46, 8}, {46, -70}, {-142, -70}, {-142, -8}, {-124, -8}}, color = {0, 0, 127}));
    connect(Tset.y, compCtrl.T_set) annotation(
      Line(points = {{-134, 4}, {-124, 4}}, color = {0, 0, 127}));
    connect(pipe2.port_b, BufferTank.port_a) annotation(
      Line(points = {{-14, 96}, {42, 96}, {42, 82}}, color = {0, 127, 255}));
    connect(pipe1.port_b, jun.port_1) annotation(
      Line(points = {{-74, 56}, {-74, 62}}, color = {0, 127, 255}));
    connect(jun.port_2, pipe2.port_a) annotation(
      Line(points = {{-74, 82}, {-74, 96}, {-34, 96}}, color = {0, 127, 255}));
    connect(jun.port_3, val.port_a) annotation(
      Line(points = {{-64, 72}, {-54, 72}}, color = {0, 127, 255}));
    connect(TsetAir.y, pumpCtrl.T_set) annotation(
      Line(points = {{101, 124}, {75, 124}, {75, 114}}, color = {0, 0, 127}));
    connect(roomWithoutLossnay1.T, pumpCtrl.T_Air) annotation(
      Line(points = {{137, 72}, {150, 72}, {150, 134}, {64, 134}, {64, 116}, {66, 116}, {66, 114}, {64, 114}}, color = {0, 0, 127}));
    connect(roomCycle1ZoneTempCtrl1.port_con, roomWithoutLossnay1.HeatPortCon) annotation(
      Line(points = {{100, 78}, {124, 78}}, color = {191, 0, 0}));
    connect(roomCycle1ZoneTempCtrl1.port_rad, roomWithoutLossnay1.HeatPortRad) annotation(
      Line(points = {{100, 68}, {124, 68}}, color = {191, 0, 0}));
    connect(pumpCtrl.y, roomCycle1ZoneTempCtrl1.m_flow_rate) annotation(
      Line(points = {{70, 94}, {70, 82}, {78, 82}}, color = {0, 0, 127}));
    connect(BufferTank.fluPorVol[1], roomCycle1ZoneTempCtrl1.port_a) annotation(
      Line(points = {{30, 82}, {66, 82}, {66, 78}, {80, 78}}, color = {0, 127, 255}));
    connect(BufferTank.fluPorVol[2], roomCycle1ZoneTempCtrl1.port_b) annotation(
      Line(points = {{30, 82}, {66, 82}, {66, 68}, {80, 68}}, color = {0, 127, 255}));
    annotation(
      Diagram(coordinateSystem(extent = {{-160, 160}, {200, -100}})),
      Documentation(info = "<html><p>
Air To Waterのうち，水回路と部屋をあわせたモデル．
</p></html>", revisions = "<html>
<ul>
<li>
</li>
<li>
</li>
<li>
</li>
</ul>
</html>"),
      experiment(StartTime = 0, StopTime = 14400, Tolerance = 1e-06, Interval = 28.8));
  end plant;
  annotation(
    uses(Modelica(version = "4.0.0"), Buildings(version = "10.0.0")));
end System;