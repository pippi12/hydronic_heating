package System
  package MyComponents
    model HeatPump
      import SI = Modelica.Units.SI;
      replaceable package Medium = Buildings.Media.Water "Medium in the pipe";
      inner Modelica.Fluid.System system annotation(
        Placement(visible = true, transformation(origin = {-128, 88}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      // parameters
      parameter Real K = 150 "1st delay K";
      parameter Modelica.Units.SI.Time T = 500 "1st delay T";
      parameter Modelica.Units.SI.Time L = 100 "wasted time L";
      parameter Modelica.Units.SI.Temperature InitialTemp = 273.15 + 30 "Initial Temperature" annotation(
        Dialog(tab = "Initialization"));
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

model TwoPipesCounterFlow
  replaceable package Medium1 = Buildings.Media.Water;
  replaceable package Medium2 = Buildings.Media.Water;
//  final redeclare Medium1 = Medium;
  extends Buildings.Fluid.Interfaces.PartialFourPortInterface(
    final show_T = false);
  parameter Modelica.Units.SI.Density rho1_default=Medium1.density_pTX(
      p=Medium1.p_default,
      T=Medium1.T_default,
      X=Medium1.X_default)
    "Default density of supply side (e.g., rho_liquidWater = 995, rho_air = 1.2)"
    annotation (Dialog(group="Advanced"));
  parameter Modelica.Units.SI.Density rho2_default=Medium2.density_pTX(
      p=Medium2.p_default,
      T=Medium2.T_default,
      X=Medium2.X_default)
    "Default density of return side (e.g., rho_liquidWater = 995, rho_air = 1.2)"
    annotation (Dialog(group="Advanced"));

  //parameter Modelica.Units.SI.Temperature Amb_T = 273.15 + 5 "Ambience temperature";
  parameter Boolean from_dp = false "= true, use m_flow = f(dp) else dp = f(m_flow)" annotation(
    Dialog(tab = "Advanced"));
  parameter Boolean have_pipCap = true "= true, a mixing volume is added that corresponds
      to the heat capacity of the pipe wall" annotation(
    Dialog(tab = "Advanced"));
  parameter Boolean have_symmetry = true "= false, the mixing volume is only on port_b,
      which improve performances, but reduces dynamic accuracy." annotation(
    Dialog(tab = "Advanced"));
  parameter Real ReC = 4000 "Reynolds number where transition to turbulence starts";
  parameter Real fac = 1 "Factor to take into account flow resistance of bends etc.,
      fac=dp_nominal/dpStraightPipe_nominal";
  parameter Modelica.Units.SI.Velocity v_nominal = 1.5 "Velocity at m_flow_nominal (used to compute default value for hydraulic diameter dh)" annotation(
    Dialog(group = "Nominal condition"));
  parameter Modelica.Units.SI.Length dh1 = sqrt(4 * m1_flow_nominal / rho1_default / v_nominal / Modelica.Constants.pi) "Hydraulic diameter of Supply pipe (assuming a round cross section area)" annotation(
    Dialog(group = "Material"));
  parameter Modelica.Units.SI.Length dh2 = sqrt(4 * m2_flow_nominal / rho2_default / v_nominal / Modelica.Constants.pi) "Hydraulic diameter of Return pipe (assuming a round cross section area)" annotation(
    Dialog(group = "Material"));
  parameter Modelica.Units.SI.Height roughness = 2.5e-5 "Average height of surface asperities (default: smooth steel pipe)" annotation(
    Dialog(group = "Material"));
  parameter Modelica.Units.SI.Length length "Pipe length" annotation(
    Dialog(group = "Material"));
  parameter Modelica.Units.SI.SpecificHeatCapacity cPip = 2300 "Specific heat of pipe wall material. 2300 for PE, 500 for steel" annotation(
    Dialog(group = "Material"));
  //, enable=have_pipCap));
  parameter Modelica.Units.SI.Density rhoPip(displayUnit = "kg/m3") = 930 "Density of pipe wall material. 930 for PE, 8000 for steel" annotation(
    Dialog(group = "Material"));
  //, enable=have_pipCap));
  parameter Modelica.Units.SI.Length thickness = 0.0035 "Pipe wall thickness" annotation(
    Dialog(group = "Material"));
  parameter Modelica.Units.SI.Length dIns "Thickness of pipe insulation, used to compute R" annotation(
    Dialog(group = "Thermal resistance"));
  parameter Modelica.Units.SI.ThermalConductivity kIns "Heat conductivity of pipe insulation, used to compute R" annotation(
    Dialog(group = "Thermal resistance"));
  parameter Modelica.Units.SI.Temperature T1_start_in(start = Medium1.T_default) = Medium1.T_default "Initialization temperature at supply pipe inlet" annotation(
    Dialog(tab = "Initialization"));
  parameter Modelica.Units.SI.Temperature T1_start_out(start = Medium1.T_default) = T1_start_in "Initialization temperature at supply pipe outlet" annotation(
    Dialog(tab = "Initialization"));
  parameter Modelica.Units.SI.Temperature T2_start_in(start = Medium2.T_default) = Medium2.T_default "Initialization temperature at return pipe inlet" annotation(
    Dialog(tab = "Initialization"));
  parameter Modelica.Units.SI.Temperature T2_start_out(start = Medium2.T_default) = T2_start_in "Initialization temperature at return pipe outlet" annotation(
    Dialog(tab = "Initialization"));
  parameter Boolean initDelay = false "Initialize delay for a constant mass flow rate if true, otherwise start from 0" annotation(
    Dialog(tab = "Initialization"));
  parameter Modelica.Units.SI.MassFlowRate m_flow_start = 0 "Initial value of mass flow rate through pipe" annotation(
    Dialog(tab = "Initialization"));
  //, enable=initDelay));
  parameter Real R1(unit = "(m.K)/W") = 1 / (kIns * 2 * Modelica.Constants.pi / Modelica.Math.log((dh1 / 2 + thickness + dIns) / (dh1 / 2 + thickness))) "Thermal resistance per unit length from fluid to boundary temperature (Supply pipe)" annotation(
    Dialog(group = "Thermal resistance"));
  parameter Real R2(unit = "(m.K)/W") = 1 / (kIns * 2 * Modelica.Constants.pi / Modelica.Math.log((dh2 / 2 + thickness + dIns) / (dh2 / 2 + thickness))) "Thermal resistance per unit length from fluid to boundary temperature (Return pipe)" annotation(
    Dialog(group = "Thermal resistance"));
  parameter Boolean linearized = false "= true, use linear relation between m_flow and dp for any flow rate" annotation(
    Evaluate = true,
    Dialog(tab = "Advanced"));
  //
  //Components
  Buildings.Fluid.FixedResistances.BaseClasses.PlugFlowPipe pipSup(
    redeclare package Medium = Medium1,
    R = R1, ReC = ReC, T_start_in = T1_start_in, T_start_out = T1_start_out, allowFlowReversal = allowFlowReversal1, cPip = cPip, dIns = dIns, dh = dh1, fac = fac, from_dp = from_dp, have_pipCap = have_pipCap, have_symmetry = have_symmetry, initDelay = initDelay, kIns = kIns, length = length, linearized = linearized, m_flow_nominal = m1_flow_nominal, m_flow_small = m1_flow_small, m_flow_start = m_flow_start, rhoPip = rhoPip, roughness = roughness, show_T = false, thickness = thickness, v_nominal = v_nominal) annotation(
    Placement(visible = true, transformation(origin = {0, 60}, extent = {{-10, 10}, {10, -10}}, rotation = 0)));
  Buildings.Fluid.FixedResistances.BaseClasses.PlugFlowPipe pipRet(
    redeclare package Medium = Medium2,
    R = R2, ReC = ReC, T_start_in = T2_start_in, T_start_out = T2_start_out, allowFlowReversal = allowFlowReversal2, cPip = cPip, dIns = dIns, dh = dh2, fac = fac, from_dp = from_dp, have_pipCap = have_pipCap, have_symmetry = have_symmetry, initDelay = initDelay, kIns = kIns, length = length, linearized = linearized, m_flow_nominal = m2_flow_nominal, m_flow_small = m2_flow_small, m_flow_start = m_flow_start, rhoPip = rhoPip, roughness = roughness, show_T = false, thickness = thickness, v_nominal = v_nominal) annotation(
    Placement(visible = true, transformation(origin = {0, -60}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPort annotation(
        Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {1, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));

equation
      connect(port_a2, pipRet.port_a) annotation(
        Line(points = {{100, -60}, {10, -60}}));
      connect(pipRet.port_b, port_b2) annotation(
        Line(points = {{-10, -60}, {-100, -60}}, color = {0, 127, 255}));
      connect(pipSup.port_b, port_b1) annotation(
        Line(points = {{10, 60}, {100, 60}}, color = {0, 127, 255}));
      connect(pipSup.port_a, port_a1) annotation(
        Line(points = {{-10, 60}, {-100, 60}}, color = {0, 127, 255}));
  connect(heatPort, pipSup.heatPort) annotation(
        Line(points = {{0, 0}, {0, 50}}, color = {191, 0, 0}));
  connect(heatPort, pipRet.heatPort) annotation(
        Line(points = {{0, 0}, {0, -50}}, color = {191, 0, 0}));
      annotation(
    Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-100, -100}, {100, 100}}), graphics = {Rectangle(origin = {0, 60.5}, fillColor = {192, 192, 192}, fillPattern = FillPattern.HorizontalCylinder, extent = {{-100, 31}, {100, -30}}), Rectangle(origin = {0, 59.5}, fillColor = {255, 0, 0}, fillPattern = FillPattern.HorizontalCylinder, extent = {{-100, 21}, {100, -20}}), Rectangle(origin = {0, 50}, lineColor = {175, 175, 175}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Backward, extent = {{-100, 50}, {100, 40}}), Rectangle(origin = {0, 70}, lineColor = {175, 175, 175}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Backward, extent = {{-100, -40}, {100, -50}}), Rectangle(origin = {0, 60}, fillColor = {215, 202, 187}, fillPattern = FillPattern.HorizontalCylinder, extent = {{-30, 20}, {28, -20}}), Text(origin = {3, -27}, extent = {{-102, -76}, {98, -104}}, textString = "L = %length"), Text(origin = {0, -17}, extent = {{-100, -56}, {100, -74}}, textString = "L = %length"), Rectangle(origin = {0, -70}, lineColor = {175, 175, 175}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Backward, extent = {{-100, 50}, {100, 40}}), Rectangle(origin = {0, -60.5}, fillColor = {192, 192, 192}, fillPattern = FillPattern.HorizontalCylinder, extent = {{-100, 31}, {100, -30}}), Rectangle(origin = {0, -59.5}, fillColor = {0, 127, 255}, fillPattern = FillPattern.HorizontalCylinder, extent = {{-100, 21}, {100, -20}}), Rectangle(origin = {0, -59}, fillColor = {215, 202, 187}, fillPattern = FillPattern.HorizontalCylinder, extent = {{-30, 20}, {28, -20}}), Rectangle(origin = {0, -50}, lineColor = {175, 175, 175}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Backward, extent = {{-100, -40}, {100, -50}})}));

end TwoPipesCounterFlow;
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

    model twoPipesCounterFlow
      extends Modelica.Icons.Example;
      replaceable package Medium = Buildings.Media.Water "Medium in the pipe" annotation(
        choicesAllMatching = true);
      final parameter Modelica.Units.SI.MassFlowRate m_flow_nominal = 6 "Mass flow rate";
      Modelica.Blocks.Sources.Ramp Tin(height = 20, duration = 0, offset = 273.15 + 50, startTime = 100) "Ramp pressure signal" annotation(
        Placement(transformation(extent = {{-100, -10}, {-80, 10}})));
      Buildings.Fluid.Sources.Boundary_pT sin(redeclare package Medium = Medium, T = 273.15 + 10, nPorts = 2, p(displayUnit = "Pa") = 101325) "Pressure boundary condition" annotation(
        Placement(visible = true, transformation(extent = {{98, 34}, {78, 54}}, rotation = 0)));
      Buildings.Fluid.Sources.MassFlowSource_T sou(redeclare package Medium = Medium, m_flow = 0.01, nPorts = 2, use_T_in = true) "Flow source" annotation(
        Placement(visible = true, transformation(extent = {{-62, 40}, {-42, 60}}, rotation = 0)));
      Buildings.Fluid.Sensors.TemperatureTwoPort senTemOut(redeclare package Medium = Medium, m_flow_nominal = m_flow_nominal, tau = 0, T_start = 323.15) "Temperature sensor" annotation(
        Placement(transformation(extent = {{40, 10}, {60, 30}})));
      Buildings.Fluid.Sensors.TemperatureTwoPort senTemIn(redeclare package Medium = Medium, m_flow_nominal = m_flow_nominal, tau = 0, T_start = 323.15) "Temperature sensor" annotation(
        Placement(transformation(extent = {{-30, 10}, {-10, 30}})));
      inner Modelica.Fluid.System system annotation(
        Placement(transformation(origin = {-90, 68}, extent = {{-10, -10}, {10, 10}})));
      System.MyComponents.TwoPipesCounterFlow twoPipesCounterFlow(
        redeclare package Medium1 = Medium,
        redeclare package Medium2 = Medium,
        allowFlowReversal1 = false, allowFlowReversal2 = false, dIns = 1, kIns = 0.001, length = 10, m1_flow_nominal = 0.1, m2_flow_nominal = 0.2) annotation(
        Placement(visible = true, transformation(origin = {12, 14}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Buildings.HeatTransfer.Sources.FixedTemperature outdoor(T = 273.15 + 5) "Boundary temperature" annotation(
        Placement(visible = true, transformation(origin = {-50, -5}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
    
    equation
  connect(Tin.y, sou.T_in) annotation(
        Line(points = {{-79, 0}, {-68, 0}, {-68, 54}, {-64, 54}}, color = {0, 0, 127}));
      connect(senTemOut.port_b, sin.ports[1]) annotation(
        Line(points = {{60, 20}, {70, 20}, {70, 44}, {78, 44}}, color = {0, 127, 255}));
  connect(sou.ports[1], senTemIn.port_a) annotation(
        Line(points = {{-42, 50}, {-35, 50}, {-35, 20}, {-30, 20}}, color = {0, 127, 255}));
  connect(senTemIn.port_b, twoPipesCounterFlow.port_a1) annotation(
        Line(points = {{-10, 20}, {2, 20}}, color = {0, 127, 255}));
  connect(twoPipesCounterFlow.port_b1, senTemOut.port_a) annotation(
        Line(points = {{22, 20}, {40, 20}}, color = {0, 127, 255}));
  connect(sou.ports[2], twoPipesCounterFlow.port_a2) annotation(
        Line(points = {{-42, 50}, {30, 50}, {30, 8}, {22, 8}}, color = {0, 127, 255}));
  connect(twoPipesCounterFlow.port_b2, sin.ports[2]) annotation(
        Line(points = {{2, 8}, {-4, 8}, {-4, -6}, {74, -6}, {74, 44}, {78, 44}}, color = {0, 127, 255}));
  connect(outdoor.port, twoPipesCounterFlow.heatPort) annotation(
        Line(points = {{-42, -4}, {-8, -4}, {-8, 14}, {12, 14}}, color = {191, 0, 0}));
      annotation(
        experiment(StartTime = 0, StopTime = 1000, Tolerance = 1e-6, Interval = 2),
        Diagram(coordinateSystem(extent = {{-100, 80}, {100, -20}})));
    end twoPipesCounterFlow;
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

  model plantEx
    extends Modelica.Icons.Example;
    replaceable package Medium = Buildings.Media.Water;
    inner Modelica.Fluid.System system annotation(
      Placement(visible = true, transformation(origin = {-88, 88}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    // Parameters
    parameter Modelica.Units.SI.Temperature T_amb = 273.15 + 5 "Ambience temperature";
    // Initialization
    parameter Modelica.Units.SI.Temperature T_medium_start = 273.15 + 30 "Water initial temperature" annotation(
      Dialog(tab = "Initialization"));
    parameter Modelica.Units.SI.Power q_flow_nominal = 5000 "Rated heat dissipation amount" annotation(
      Dialog(tab = "Initialization"));
    parameter Modelica.Units.SI.Temperature T_rad_start = 273.15 + 40 "Radiator water initial temperature" annotation(
      Dialog(tab = "Initialization"));
    // Pipe
    //parameter Modelica.Units.SI.Length pip_len = 1 "Length of pipe";
    parameter Modelica.Units.SI.SpecificHeatCapacity pip_c = 386 "Specific heat of pipe wall material" annotation(
      Dialog(group = "Pipe"));
    parameter Modelica.Units.SI.Density pip_rho = 8960 "Density of pipe wall material" annotation(
      Dialog(group = "Pipe"));
    parameter Modelica.Units.SI.Length pip_dh = 0.025 "Inner diameter of pipe" annotation(
      Dialog(group = "Pipe"));
    parameter Modelica.Units.SI.Length pip_thickness = 0.005 "Thickness of pipe wall" annotation(
      Dialog(group = "Pipe"));
    parameter Modelica.Units.SI.ThermalConductivity pip_kIns = 398 "Heat conductivity of pipe insulation" annotation(
      Dialog(group = "Pipe"));
    parameter Modelica.Units.SI.CoefficientOfHeatTransfer pip_air_htc = 10 "Air-Pipe heat transfer coefficient" annotation(
      Dialog(group = "Pipe"));
    parameter Modelica.Units.SI.MassFlowRate m_flow_nominal = 0.1 "Mass flow rate" annotation(
      Dialog(group = "Pipe"));
    // Pump
    parameter Real m_flow = 10 / 60 * 0.001 * 1000 "Mass flow" annotation(
      Dialog(group = "Pump"));
    // ===
    // Components
    // ===
    Buildings.HeatTransfer.Sources.FixedTemperature outdoor(T = T_amb) "Boundary temperature" annotation(
      Placement(visible = true, transformation(origin = {-110, 45}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
    Buildings.Fluid.Movers.FlowControlled_m_flow pump_m_flow(
      redeclare package Medium = Medium, T_start = T_medium_start, allowFlowReversal = false, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, inputType = Buildings.Fluid.Types.InputType.Continuous, m_flow_nominal = m_flow_nominal, massFlowRates = {0, 0.5, 1} * m_flow_nominal, nominalValuesDefineDefaultPressureCurve = true, use_inputFilter = false) "Pump with m_flow input" annotation(
      Placement(visible = true, transformation(origin = {-54, 62}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Ramp ramp1(duration = 0, height = m_flow, offset = 0, startTime = 0) annotation(
      Placement(visible = true, transformation(origin = {-35, 89}, extent = {{7, -7}, {-7, 7}}, rotation = 0)));
    System.MyComponents.HeatPump heatPump(
      redeclare package Medium = Medium, InitialTemp = T_medium_start) annotation(
      Placement(visible = true, transformation(origin = {-74, -8}, extent = {{-20, 18}, {10, 44}}, rotation = 90)));
    System.MyComponents.TwoPipesCounterFlow pips1(
      redeclare package Medium1 = Medium,
      redeclare package Medium2 = Medium, R1 = Modelica.Math.log((pip_dh + pip_thickness) / pip_dh) / (2 * Modelica.Constants.pi * pip_kIns) + 1 / pip_air_htc / Modelica.Constants.pi / (pip_dh + pip_thickness), R2 = Modelica.Math.log((pip_dh + pip_thickness) / pip_dh) / (2 * Modelica.Constants.pi * pip_kIns) + 1 / pip_air_htc / Modelica.Constants.pi / (pip_dh + pip_thickness), T1_start_in = T_medium_start, T1_start_out = T_medium_start, T2_start_in = T_medium_start, T2_start_out = T_medium_start, allowFlowReversal1 = false, allowFlowReversal2 = false, cPip = pip_c, dIns = pip_thickness, dh1 = pip_dh, dh2 = pip_dh, kIns = pip_kIns, length = 1, m1_flow_nominal = m_flow_nominal, m2_flow_nominal = m_flow_nominal, m_flow_start = 0, rhoPip = pip_rho, thickness = pip_thickness) annotation(
        Placement(visible = true, transformation(origin = {-68, 24}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
    System.MyComponents.TwoPipesCounterFlow pips_kit_1(
      redeclare package Medium1 = Medium,
      redeclare package Medium2 = Medium, R1 = Modelica.Math.log((pip_dh + pip_thickness) / pip_dh) / (2 * Modelica.Constants.pi * pip_kIns) + 1 / pip_air_htc / Modelica.Constants.pi / (pip_dh + pip_thickness), R2 = Modelica.Math.log((pip_dh + pip_thickness) / pip_dh) / (2 * Modelica.Constants.pi * pip_kIns) + 1 / pip_air_htc / Modelica.Constants.pi / (pip_dh + pip_thickness), T1_start_in = T_medium_start, T1_start_out = T_medium_start, T2_start_in = T_medium_start, T2_start_out = T_medium_start, allowFlowReversal1 = false, allowFlowReversal2 = false, cPip = pip_c, dIns = pip_thickness, dh1 = pip_dh, dh2 = pip_dh, kIns = pip_kIns, length = 1, m1_flow_nominal = m_flow_nominal, m2_flow_nominal = m_flow_nominal, m_flow_start = 0, rhoPip = pip_rho, thickness = pip_thickness) annotation(
        Placement(visible = true, transformation(origin = {-40, 24}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Buildings.Fluid.FixedResistances.Junction jun_sup_kit1(
      redeclare package Medium = Medium, dp_nominal = {0, 0, 0}, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, m_flow_nominal = m_flow * {1, -1, -1}) "Splitter" annotation(
        Placement(visible = true, transformation(origin = {0, -12}, extent = {{10, -10}, {-10, 10}}, rotation = 180)));
    Buildings.Fluid.FixedResistances.Junction jun_sup_kit2(
      redeclare package Medium = Medium, dp_nominal = {0, 0, 0}, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, m_flow_nominal = m_flow * {1, -1, -1}) annotation(
        Placement(visible = true, transformation(origin = {30, -12}, extent = {{10, 10}, {-10, -10}}, rotation = 180)));
    Buildings.Fluid.FixedResistances.Junction jun_ret_kit1(
      redeclare package Medium = Medium, dp_nominal = {0, 0, 0}, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, m_flow_nominal = m_flow * {1, -1, 1}) annotation(
        Placement(visible = true, transformation(origin = {12, -34}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
    Buildings.Fluid.FixedResistances.Junction jun_ret_kit2(
      redeclare package Medium = Medium, dp_nominal = {0, 0, 0}, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, m_flow_nominal = m_flow * {1, -1, 1}) annotation(
        Placement(visible = true, transformation(origin = {42, -34}, extent = {{-10, 10}, {10, -10}}, rotation = 180)));
    Buildings.Fluid.FixedResistances.Junction jun_sup_kit3(
      redeclare package Medium = Medium, dp_nominal = {0, 0, 0}, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, m_flow_nominal = m_flow * {1, -1, -1}) annotation(
        Placement(visible = true, transformation(origin = {0, 52}, extent = {{10, 10}, {-10, -10}}, rotation = 270)));
    Buildings.Fluid.HeatExchangers.Radiators.RadiatorEN442_2 rad_kit1(
      redeclare package Medium = Medium, Q_flow_nominal = q_flow_nominal, TAir_nominal = 293.15, T_a_nominal = 353.15, T_b_nominal = 333.15, T_start = T_rad_start, allowFlowReversal = false, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, fraRad = 0, nEle = 1, p_start = 101325) annotation(
        Placement(visible = true, transformation(origin = {70, 40}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    Buildings.Fluid.FixedResistances.Junction jun_ret_kit3(
      redeclare package Medium = Medium, dp_nominal = {0, 0, 0}, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, m_flow_nominal = m_flow * {1, -1, 1}) annotation(
        Placement(visible = true, transformation(origin = {12, 40}, extent = {{-10, 10}, {10, -10}}, rotation = 270)));
    System.MyComponents.TwoPipesCounterFlow pips_kit_2(
      redeclare package Medium1 = Medium,
      redeclare package Medium2 = Medium, R1 = Modelica.Math.log((pip_dh + pip_thickness) / pip_dh) / (2 * Modelica.Constants.pi * pip_kIns) + 1 / pip_air_htc / Modelica.Constants.pi / (pip_dh + pip_thickness), R2 = Modelica.Math.log((pip_dh + pip_thickness) / pip_dh) / (2 * Modelica.Constants.pi * pip_kIns) + 1 / pip_air_htc / Modelica.Constants.pi / (pip_dh + pip_thickness), T1_start_in = T_medium_start, T1_start_out = T_medium_start, T2_start_in = T_medium_start, T2_start_out = T_medium_start, allowFlowReversal1 = false, allowFlowReversal2 = false, cPip = pip_c, dIns = pip_thickness, dh1 = pip_dh, dh2 = pip_dh, kIns = pip_kIns, length = 1, m1_flow_nominal = m_flow_nominal, m2_flow_nominal = m_flow_nominal, m_flow_start = 0, rhoPip = pip_rho, thickness = pip_thickness) annotation(
        Placement(visible = true, transformation(origin = {6, 14}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
    System.MyComponents.TwoPipesCounterFlow pips_kit_3(
      redeclare package Medium1 = Medium,
      redeclare package Medium2 = Medium, R1 = Modelica.Math.log((pip_dh + pip_thickness) / pip_dh) / (2 * Modelica.Constants.pi * pip_kIns) + 1 / pip_air_htc / Modelica.Constants.pi / (pip_dh + pip_thickness), R2 = Modelica.Math.log((pip_dh + pip_thickness) / pip_dh) / (2 * Modelica.Constants.pi * pip_kIns) + 1 / pip_air_htc / Modelica.Constants.pi / (pip_dh + pip_thickness), T1_start_in = T_medium_start, T1_start_out = T_medium_start, T2_start_in = T_medium_start, T2_start_out = T_medium_start, allowFlowReversal1 = false, allowFlowReversal2 = false, cPip = pip_c, dIns = pip_thickness, dh1 = pip_dh, dh2 = pip_dh, kIns = pip_kIns, length = 1, m1_flow_nominal = m_flow_nominal, m2_flow_nominal = m_flow_nominal, m_flow_start = 0.01, rhoPip = pip_rho, thickness = pip_thickness) annotation(
        Placement(visible = true, transformation(origin = {42, 46}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    System.MyComponents.TwoPipesCounterFlow pips_kit_4(
      redeclare package Medium1 = Medium,
      redeclare package Medium2 = Medium, R1 = Modelica.Math.log((pip_dh + pip_thickness) / pip_dh) / (2 * Modelica.Constants.pi * pip_kIns) + 1 / pip_air_htc / Modelica.Constants.pi / (pip_dh + pip_thickness), R2 = Modelica.Math.log((pip_dh + pip_thickness) / pip_dh) / (2 * Modelica.Constants.pi * pip_kIns) + 1 / pip_air_htc / Modelica.Constants.pi / (pip_dh + pip_thickness), T1_start_in = T_medium_start, T1_start_out = T_medium_start, T2_start_in = T_medium_start, T2_start_out = T_medium_start, allowFlowReversal1 = false, allowFlowReversal2 = false, cPip = pip_c, dIns = pip_thickness, dh1 = pip_dh, dh2 = pip_dh, kIns = pip_kIns, length = 1, m1_flow_nominal = m_flow_nominal, m2_flow_nominal = m_flow_nominal, m_flow_start = 0, rhoPip = pip_rho, thickness = pip_thickness) annotation(
        Placement(visible = true, transformation(origin = {40, 88}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Buildings.Fluid.HeatExchangers.Radiators.RadiatorEN442_2 rad_liv1(
      redeclare package Medium = Medium, Q_flow_nominal = q_flow_nominal, TAir_nominal = 293.15, T_a_nominal = 353.15, T_b_nominal = 333.15, T_start = T_rad_start, allowFlowReversal = false, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, fraRad = 0, nEle = 1, p_start = 101325) annotation(
        Placement(visible = true, transformation(origin = {98, 88}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    System.MyComponents.TwoPipesCounterFlow pips_hal_1(
      redeclare package Medium1 = Medium,
      redeclare package Medium2 = Medium, R1 = Modelica.Math.log((pip_dh + pip_thickness) / pip_dh) / (2 * Modelica.Constants.pi * pip_kIns) + 1 / pip_air_htc / Modelica.Constants.pi / (pip_dh + pip_thickness), R2 = Modelica.Math.log((pip_dh + pip_thickness) / pip_dh) / (2 * Modelica.Constants.pi * pip_kIns) + 1 / pip_air_htc / Modelica.Constants.pi / (pip_dh + pip_thickness), T1_start_in = T_medium_start, T1_start_out = T_medium_start, T2_start_in = T_medium_start, T2_start_out = T_medium_start, allowFlowReversal1 = false, allowFlowReversal2 = false, cPip = pip_c, dIns = pip_thickness, dh1 = pip_dh, dh2 = pip_dh, kIns = pip_kIns, length = 1, m1_flow_nominal = m_flow_nominal, m2_flow_nominal = m_flow_nominal, m_flow_start = 0, rhoPip = pip_rho, thickness = pip_thickness) annotation(
        Placement(visible = true, transformation(origin = {36, -62}, extent = {{10, -10}, {-10, 10}}, rotation = 90)));
    Buildings.Fluid.FixedResistances.Junction jun_ret_hal1(
      redeclare package Medium = Medium, dp_nominal = {0, 0, 0}, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, m_flow_nominal = m_flow * {1, -1, 1}) annotation(
        Placement(visible = true, transformation(origin = {42, -110}, extent = {{10, 10}, {-10, -10}}, rotation = 270)));
    Buildings.Fluid.FixedResistances.Junction jun_sup_hal1(
      redeclare package Medium = Medium, dp_nominal = {0, 0, 0}, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, m_flow_nominal = m_flow * {1, -1, -1}) annotation(
        Placement(visible = true, transformation(origin = {30, -98}, extent = {{-10, 10}, {10, -10}}, rotation = 270)));
    Buildings.Fluid.HeatExchangers.Radiators.RadiatorEN442_2 rad_hal1(
      redeclare package Medium = Medium, Q_flow_nominal = q_flow_nominal, TAir_nominal = 293.15, T_a_nominal = 353.15, T_b_nominal = 333.15, T_start = T_rad_start, allowFlowReversal = false, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, fraRad = 0, nEle = 1, p_start = 101325) annotation(
        Placement(visible = true, transformation(origin = {112, -112}, extent = {{10, -10}, {-10, 10}}, rotation = 90)));
    Buildings.Fluid.HeatExchangers.Radiators.RadiatorEN442_2 rad_toi(
      redeclare package Medium = Medium, Q_flow_nominal = q_flow_nominal, TAir_nominal = 293.15, T_a_nominal = 353.15, T_b_nominal = 333.15, T_start = T_rad_start, allowFlowReversal = false, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, fraRad = 0, nEle = 1, p_start = 101325) annotation(
        Placement(visible = true, transformation(origin = {-50, -164}, extent = {{10, -10}, {-10, 10}}, rotation = -90)));
    System.MyComponents.TwoPipesCounterFlow pips_hal_2(
      redeclare package Medium1 = Medium,
      redeclare package Medium2 = Medium, R1 = Modelica.Math.log((pip_dh + pip_thickness) / pip_dh) / (2 * Modelica.Constants.pi * pip_kIns) + 1 / pip_air_htc / Modelica.Constants.pi / (pip_dh + pip_thickness), R2 = Modelica.Math.log((pip_dh + pip_thickness) / pip_dh) / (2 * Modelica.Constants.pi * pip_kIns) + 1 / pip_air_htc / Modelica.Constants.pi / (pip_dh + pip_thickness), T1_start_in = T_medium_start, T1_start_out = T_medium_start, T2_start_in = T_medium_start, T2_start_out = T_medium_start, allowFlowReversal1 = false, allowFlowReversal2 = false, cPip = pip_c, dIns = pip_thickness, dh1 = pip_dh, dh2 = pip_dh, kIns = pip_kIns, length = 1, m1_flow_nominal = m_flow_nominal, m2_flow_nominal = m_flow_nominal, m_flow_start = 0, rhoPip = pip_rho, thickness = pip_thickness) annotation(
        Placement(visible = true, transformation(origin = {4, -184}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    System.MyComponents.TwoPipesCounterFlow pips_toi_1(
      redeclare package Medium1 = Medium,
      redeclare package Medium2 = Medium, R1 = Modelica.Math.log((pip_dh + pip_thickness) / pip_dh) / (2 * Modelica.Constants.pi * pip_kIns) + 1 / pip_air_htc / Modelica.Constants.pi / (pip_dh + pip_thickness), R2 = Modelica.Math.log((pip_dh + pip_thickness) / pip_dh) / (2 * Modelica.Constants.pi * pip_kIns) + 1 / pip_air_htc / Modelica.Constants.pi / (pip_dh + pip_thickness), T1_start_in = T_medium_start, T1_start_out = T_medium_start, T2_start_in = T_medium_start, T2_start_out = T_medium_start, allowFlowReversal1 = false, allowFlowReversal2 = false, cPip = pip_c, dIns = pip_thickness, dh1 = pip_dh, dh2 = pip_dh, kIns = pip_kIns, length = 1, m1_flow_nominal = m_flow_nominal, m2_flow_nominal = m_flow_nominal, m_flow_start = 0, rhoPip = pip_rho, thickness = pip_thickness) annotation(
        Placement(visible = true, transformation(origin = {-24, -184}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    Buildings.Fluid.FixedResistances.Junction jun_sup_hal2(
      redeclare package Medium = Medium, dp_nominal = {0, 0, 0}, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, m_flow_nominal = m_flow * {1, -1, -1}) annotation(
        Placement(visible = true, transformation(origin = {30, -158}, extent = {{-10, 10}, {10, -10}}, rotation = 270)));
    Buildings.Fluid.FixedResistances.Junction jun_ret_hal2(
      redeclare package Medium = Medium, dp_nominal = {0, 0, 0}, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, m_flow_nominal = m_flow * {1, -1, 1}) annotation(
        Placement(visible = true, transformation(origin = {42, -146}, extent = {{10, 10}, {-10, -10}}, rotation = 270)));
    Buildings.Fluid.HeatExchangers.Radiators.RadiatorEN442_2 rad_liv2(
      redeclare package Medium = Medium, Q_flow_nominal = q_flow_nominal, TAir_nominal = 293.15, T_a_nominal = 353.15, T_b_nominal = 333.15, T_start = T_rad_start, allowFlowReversal = false, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, fraRad = 0, nEle = 1, p_start = 101325) annotation(
        Placement(visible = true, transformation(origin = {150, -158}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    System.MyComponents.TwoPipesCounterFlow pips_kit_31(
      redeclare package Medium1 = Medium,
      redeclare package Medium2 = Medium, R1 = Modelica.Math.log((pip_dh + pip_thickness) / pip_dh) / (2 * Modelica.Constants.pi * pip_kIns) + 1 / pip_air_htc / Modelica.Constants.pi / (pip_dh + pip_thickness), R2 = Modelica.Math.log((pip_dh + pip_thickness) / pip_dh) / (2 * Modelica.Constants.pi * pip_kIns) + 1 / pip_air_htc / Modelica.Constants.pi / (pip_dh + pip_thickness), T1_start_in = T_medium_start, T1_start_out = T_medium_start, T2_start_in = T_medium_start, T2_start_out = T_medium_start, allowFlowReversal1 = false, allowFlowReversal2 = false, cPip = pip_c, dIns = pip_thickness, dh1 = pip_dh, dh2 = pip_dh, kIns = pip_kIns, length = 1, m1_flow_nominal = m_flow_nominal, m2_flow_nominal = m_flow_nominal, m_flow_start = 0, rhoPip = pip_rho, thickness = pip_thickness) annotation(
        Placement(visible = true, transformation(origin = {74, -104}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    System.MyComponents.TwoPipesCounterFlow pips_liv_1(
      redeclare package Medium1 = Medium,
      redeclare package Medium2 = Medium, R1 = Modelica.Math.log((pip_dh + pip_thickness) / pip_dh) / (2 * Modelica.Constants.pi * pip_kIns) + 1 / pip_air_htc / Modelica.Constants.pi / (pip_dh + pip_thickness), R2 = Modelica.Math.log((pip_dh + pip_thickness) / pip_dh) / (2 * Modelica.Constants.pi * pip_kIns) + 1 / pip_air_htc / Modelica.Constants.pi / (pip_dh + pip_thickness), T1_start_in = T_medium_start, T1_start_out = T_medium_start, T2_start_in = T_medium_start, T2_start_out = T_medium_start, allowFlowReversal1 = false, allowFlowReversal2 = false, cPip = pip_c, dIns = pip_thickness, dh1 = pip_dh, dh2 = pip_dh, kIns = pip_kIns, length = 1, m1_flow_nominal = m_flow_nominal, m2_flow_nominal = m_flow_nominal, m_flow_start = 0, rhoPip = pip_rho, thickness = pip_thickness) annotation(
      Placement(visible = true, transformation(origin = {94, -152}, extent = {{-10, 10}, {10, -10}}, rotation = 0)));
    Buildings.Fluid.Storage.ExpansionVessel exp(
      redeclare package Medium = Medium, p_start = 101325, V_start = 0.001, T_start = T_medium_start) annotation(
      Placement(visible = true, transformation(origin = {-68, -44}, extent = {{-10, -10}, {10, 10}}, rotation = 180)));
    Modelica.Blocks.Sources.RealExpression Tset(y = 20) annotation(
      Placement(visible = true, transformation(origin = {-118, -34}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Buildings.Fluid.FixedResistances.PlugFlowPipe pipe1(
      redeclare package Medium = Medium, R = Modelica.Math.log((pip_dh + pip_thickness) / pip_dh) / (2 * Modelica.Constants.pi * pip_kIns) + 1 / pip_air_htc / Modelica.Constants.pi / (pip_dh + pip_thickness), T_start_in = T_medium_start, T_start_out = T_medium_start, allowFlowReversal = false, cPip = pip_c, dIns = pip_thickness, dh = pip_dh, kIns = pip_kIns, length = 1, m_flow_nominal = m_flow_nominal, rhoPip = pip_rho, thickness = pip_thickness) "Outside PlugFlowPipe" annotation(
      Placement(visible = true, transformation(origin = {84, -24}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Buildings.Fluid.MixingVolumes.MixingVolume roomAir(redeclare package Medium = Medium, T_start = 273.15 + 10, V = 100, allowFlowReversal = true, m_flow_nominal = 0.01, massDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, nPorts = 0, p_start(displayUnit = "Pa") = 101325, use_C_flow = false) annotation(
      Placement(visible = true, transformation(origin = {208, -14}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  
  equation
    connect(heatPump.port_b, pips1.port_a1) annotation(
      Line(points = {{-74, 0}, {-74, 14}}, color = {0, 127, 255}));
    connect(pips_kit_1.port_b2, pips1.port_a2) annotation(
      Line(points = {{-46, 34}, {-46, 46}, {-62, 46}, {-62, 34}}, color = {0, 127, 255}));
    connect(pips_kit_1.port_b1, jun_sup_kit1.port_1) annotation(
      Line(points = {{-34, 14}, {-34, -12}, {-10, -12}}, color = {0, 127, 255}));
    connect(jun_sup_kit1.port_2, jun_sup_kit2.port_1) annotation(
      Line(points = {{10, -12}, {20, -12}}, color = {0, 127, 255}));
    connect(jun_ret_kit2.port_2, jun_ret_kit1.port_1) annotation(
      Line(points = {{32, -34}, {22, -34}}, color = {0, 127, 255}));
    connect(jun_ret_kit1.port_2, pips_kit_1.port_a2) annotation(
      Line(points = {{2, -34}, {-46, -34}, {-46, 14}}, color = {0, 127, 255}));
    connect(jun_sup_kit1.port_3, pips_kit_2.port_a1) annotation(
      Line(points = {{0, -2}, {0, 4}}, color = {0, 127, 255}));
    connect(pips_kit_2.port_b2, jun_ret_kit1.port_3) annotation(
      Line(points = {{12, 4}, {12, -24}}, color = {0, 127, 255}));
    connect(pips_kit_2.port_a2, jun_ret_kit3.port_2) annotation(
      Line(points = {{12, 24}, {12, 30}}, color = {0, 127, 255}));
    connect(pips_kit_2.port_b1, jun_sup_kit3.port_1) annotation(
      Line(points = {{0, 24}, {0, 42}}, color = {0, 127, 255}));
    connect(jun_sup_kit3.port_3, pips_kit_3.port_a1) annotation(
        Line(points = {{10, 52}, {32, 52}}, color = {0, 127, 255}));
    connect(pips_kit_3.port_b1, rad_kit1.port_a) annotation(
        Line(points = {{52, 52}, {80, 52}, {80, 40}}, color = {0, 127, 255}));
    connect(rad_kit1.port_b, pips_kit_3.port_a2) annotation(
        Line(points = {{60, 40}, {52, 40}}, color = {0, 127, 255}));
    connect(pips_kit_3.port_b2, jun_ret_kit3.port_3) annotation(
        Line(points = {{32, 40}, {22, 40}}, color = {0, 127, 255}));
    connect(pips_kit_4.port_b2, jun_ret_kit3.port_1) annotation(
        Line(points = {{30, 82}, {12, 82}, {12, 50}}, color = {0, 127, 255}));
    connect(jun_sup_kit3.port_2, pips_kit_4.port_a1) annotation(
        Line(points = {{0, 62}, {0, 94}, {30, 94}}, color = {0, 127, 255}));
    connect(pips_kit_4.port_b1, rad_liv1.port_a) annotation(
        Line(points = {{50, 94}, {74, 94}, {74, 98}, {98, 98}}, color = {0, 127, 255}));
    connect(rad_liv1.port_b, pips_kit_4.port_a2) annotation(
        Line(points = {{98, 78}, {74, 78}, {74, 82}, {50, 82}}, color = {0, 127, 255}));
    connect(jun_sup_kit2.port_3, pips_hal_1.port_a1) annotation(
        Line(points = {{30, -22}, {30, -52}}, color = {0, 127, 255}));
    connect(jun_ret_kit2.port_3, pips_hal_1.port_b2) annotation(
        Line(points = {{42, -44}, {42, -52}}, color = {0, 127, 255}));
    connect(jun_sup_hal1.port_3, rad_hal1.port_a) annotation(
        Line(points = {{40, -98}, {112, -98}, {112, -102}}, color = {0, 127, 255}));
    connect(jun_ret_hal1.port_2, pips_hal_1.port_a2) annotation(
        Line(points = {{42, -100}, {42, -72}}, color = {0, 127, 255}));
    connect(pips_hal_1.port_b1, jun_sup_hal1.port_1) annotation(
        Line(points = {{30, -72}, {30, -88}}, color = {0, 127, 255}));
    connect(pips_hal_2.port_b1, pips_toi_1.port_a1) annotation(
        Line(points = {{-6, -178}, {-14, -178}}, color = {0, 127, 255}));
    connect(pips_toi_1.port_b2, pips_hal_2.port_a2) annotation(
        Line(points = {{-14, -190}, {-6, -190}}, color = {0, 127, 255}));
    connect(pips_toi_1.port_b1, rad_toi.port_a) annotation(
        Line(points = {{-34, -178}, {-50, -178}, {-50, -174}}, color = {0, 127, 255}));
    connect(rad_toi.port_b, pips_toi_1.port_a2) annotation(
        Line(points = {{-50, -154}, {-50, -148}, {-60, -148}, {-60, -190}, {-34, -190}}, color = {0, 127, 255}));
    connect(jun_sup_hal1.port_2, jun_sup_hal2.port_1) annotation(
        Line(points = {{30, -108}, {30, -148}}, color = {0, 127, 255}));
    connect(jun_sup_hal2.port_2, pips_hal_2.port_a1) annotation(
        Line(points = {{30, -168}, {30, -178}, {14, -178}}, color = {0, 127, 255}));
    connect(pips_hal_2.port_b2, jun_ret_hal2.port_1) annotation(
        Line(points = {{14, -190}, {42, -190}, {42, -156}}, color = {0, 127, 255}));
    connect(jun_ret_hal2.port_2, jun_ret_hal1.port_1) annotation(
        Line(points = {{42, -136}, {42, -120}}, color = {0, 127, 255}));
    connect(jun_sup_hal2.port_3, pips_liv_1.port_a1) annotation(
        Line(points = {{40, -158}, {84, -158}}, color = {0, 127, 255}));
    connect(pips_liv_1.port_b1, rad_liv2.port_a) annotation(
        Line(points = {{104, -158}, {140, -158}}, color = {0, 127, 255}));
    connect(rad_liv2.port_b, pips_liv_1.port_a2) annotation(
        Line(points = {{160, -158}, {168, -158}, {168, -146}, {104, -146}}, color = {0, 127, 255}));
    connect(pips_liv_1.port_b2, jun_ret_hal2.port_3) annotation(
        Line(points = {{84, -146}, {52, -146}}, color = {0, 127, 255}));
  connect(pips1.port_b2, exp.port_a) annotation(
      Line(points = {{-62, 14}, {-62, -34}, {-68, -34}}, color = {0, 127, 255}));
  connect(exp.port_a, heatPump.port_a) annotation(
      Line(points = {{-68, -34}, {-74, -34}, {-74, -20}}, color = {0, 127, 255}));
    connect(pips1.port_b1, pump_m_flow.port_a) annotation(
        Line(points = {{-74, 34}, {-74, 62}, {-64, 62}}, color = {0, 127, 255}));
    connect(pump_m_flow.port_b, pips_kit_1.port_a1) annotation(
        Line(points = {{-44, 62}, {-34, 62}, {-34, 34}}, color = {0, 127, 255}));
    connect(ramp1.y, pump_m_flow.m_flow_in) annotation(
        Line(points = {{-42, 90}, {-54, 90}, {-54, 74}}, color = {0, 0, 127}));
    connect(Tset.y, heatPump.freq) annotation(
        Line(points = {{-106, -34}, {-86, -34}, {-86, -24}}, color = {0, 0, 127}));
    connect(jun_sup_kit2.port_2, pipe1.port_a) annotation(
        Line(points = {{40, -12}, {84, -12}, {84, -14}}, color = {0, 127, 255}));
    connect(pipe1.port_b, jun_ret_kit2.port_1) annotation(
      Line(points = {{84, -34}, {52, -34}}, color = {0, 127, 255}));
  connect(outdoor.port, pips1.heatPort) annotation(
      Line(points = {{-102, 46}, {-68, 46}, {-68, 24}}, color = {191, 0, 0}));
  connect(outdoor.port, pips_kit_1.heatPort) annotation(
      Line(points = {{-102, 46}, {-40, 46}, {-40, 24}}, color = {191, 0, 0}));
  connect(outdoor.port, pips_kit_4.heatPort) annotation(
      Line(points = {{-102, 46}, {40, 46}, {40, 88}}, color = {191, 0, 0}));
  connect(outdoor.port, pips_kit_3.heatPort) annotation(
      Line(points = {{-102, 46}, {42, 46}}, color = {191, 0, 0}));
  connect(outdoor.port, pips_kit_2.heatPort) annotation(
      Line(points = {{-102, 46}, {6, 46}, {6, 14}}, color = {191, 0, 0}));
  connect(outdoor.port, pipe1.heatPort) annotation(
      Line(points = {{-102, 46}, {94, 46}, {94, -24}}, color = {191, 0, 0}));
  connect(outdoor.port, pips_hal_1.heatPort) annotation(
      Line(points = {{-102, 46}, {36, 46}, {36, -62}}, color = {191, 0, 0}));
  connect(outdoor.port, pips_kit_31.heatPort) annotation(
      Line(points = {{-102, 46}, {74, 46}, {74, -104}}, color = {191, 0, 0}));
  connect(outdoor.port, pips_liv_1.heatPort) annotation(
      Line(points = {{-102, 46}, {94, 46}, {94, -152}}, color = {191, 0, 0}));
  connect(outdoor.port, pips_hal_2.heatPort) annotation(
      Line(points = {{-102, 46}, {4, 46}, {4, -184}}, color = {191, 0, 0}));
  connect(outdoor.port, pips_toi_1.heatPort) annotation(
      Line(points = {{-102, 46}, {-24, 46}, {-24, -184}}, color = {191, 0, 0}));
  connect(rad_kit1.heatPortCon, roomAir.heatPort) annotation(
      Line(points = {{72, 48}, {72, 64}, {184, 64}, {184, -14}, {198, -14}}, color = {191, 0, 0}));
  connect(rad_kit1.heatPortRad, roomAir.heatPort) annotation(
      Line(points = {{68, 48}, {68, 64}, {184, 64}, {184, -14}, {198, -14}}, color = {191, 0, 0}));
  connect(rad_liv1.heatPortCon, roomAir.heatPort) annotation(
      Line(points = {{106, 90}, {184, 90}, {184, -14}, {198, -14}}, color = {191, 0, 0}));
  connect(rad_liv1.heatPortRad, roomAir.heatPort) annotation(
      Line(points = {{106, 86}, {184, 86}, {184, -14}, {198, -14}}, color = {191, 0, 0}));
  connect(rad_hal1.heatPortCon, roomAir.heatPort) annotation(
      Line(points = {{104, -110}, {102, -110}, {102, -14}, {198, -14}}, color = {191, 0, 0}));
  connect(rad_hal1.heatPortRad, roomAir.heatPort) annotation(
      Line(points = {{104, -114}, {102, -114}, {102, -14}, {198, -14}}, color = {191, 0, 0}));
  connect(rad_liv2.heatPortCon, roomAir.heatPort) annotation(
      Line(points = {{148, -150}, {146, -150}, {146, -14}, {198, -14}}, color = {191, 0, 0}));
  connect(rad_liv2.heatPortRad, roomAir.heatPort) annotation(
      Line(points = {{152, -150}, {152, -14}, {198, -14}}, color = {191, 0, 0}));
  connect(rad_toi.heatPortRad, roomAir.heatPort) annotation(
      Line(points = {{-42, -162}, {186, -162}, {186, -14}, {198, -14}}, color = {191, 0, 0}));
  connect(rad_toi.heatPortCon, roomAir.heatPort) annotation(
      Line(points = {{-42, -166}, {186, -166}, {186, -14}, {198, -14}}, color = {191, 0, 0}));
  connect(rad_hal1.port_b, pips_kit_31.port_a2) annotation(
      Line(points = {{112, -122}, {112, -124}, {98, -124}, {98, -110}, {84, -110}}, color = {0, 127, 255}));
  connect(pips_kit_31.port_b2, jun_ret_hal1.port_3) annotation(
      Line(points = {{64, -110}, {52, -110}}, color = {0, 127, 255}));
  annotation(
      Diagram(coordinateSystem(extent = {{-120, 100}, {180, -200}})));
  end plantEx;
  annotation(
    uses(Modelica(version = "4.0.0"), Buildings(version = "10.0.0")));
end System;