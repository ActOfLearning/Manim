from manimlib.imports import *
import numpy as np
'''
Create a folder named 'assets' inside the Manim folder.
Create three folders named 'raster_images', 'sounds' and 'svg_images' inside the 'assets' folder.
'''

class Angle(VGroup):

    CONFIG = {
        'radius': 0.5,
        'color': RED,
        'opacity': 0.4,
        'stroke_width': 5,
        'below_180': True,
    }

    def __init__(self, A, O, B, **kwargs):

        VMobject.__init__(self, **kwargs)
        OA, OB = A-O, B-O
        if self.below_180:
            theta = np.angle(complex(*OA[:2])/complex(*OB[:2])) # angle of OB to OA
        else:
            theta = TAU + np.angle(complex(*OA[:2])/complex(*OB[:2]))

        if abs(abs(theta) - PI / 2) <= 1e-3:
            self.add(Elbow(angle = np.angle(complex(*OB[:2])), color = self.color, width = self.radius).shift(O))
        else:
            self.add(Sector(inner_radius = 0, outer_radius = self.radius, angle = theta, arc_center = O, color = self.color, fill_opacity = self.opacity, start_angle = Line(O, B).get_angle()))
            self.add(Arc(start_angle=Line(O, B).get_angle(), angle=theta, radius=self.radius,
                     stroke_width=self.stroke_width, color=self.color, arc_center=O))

def ParallalLine(pt0, line0):
    res = VGroup()
    if type(pt0).__name__ == "Dot" and type(line0).__name__ == "Line":
        line_start, line_end = line0.get_start_and_end()
        p_cent = pt0.get_center()
        if np.linalg.norm(line_end - line_start) > 0:
            normalized_line = (line_end - line_start) / np.linalg.norm(line_end - line_start)
            return Line(p_cent, p_cent + normalized_line)
    return res

def PerpendicularFoot(pt0, line0):
    res = VGroup()
    if type(pt0).__name__ == "Dot" and type(line0).__name__ == "Line":
        line_start, line_end = line0.get_start_and_end()
        p_cent = pt0.get_center()
        if np.linalg.norm(line_end - line_start) > 0:
            normalized_line = (line_end - line_start) / np.linalg.norm(line_end - line_start)
            line_start += np.dot(p_cent - line_start, normalized_line) * normalized_line
            return Dot().move_to(line_start)
    return res

def PerpendicularLine(pt0, line0):
    line_start, line_end = line0.get_start_and_end()
    line_unit = (line_end - line_start) / np.linalg.norm(line_end - line_start)
    perp_unit = np.cross(line_unit, OUT)
    perp_foot = PerpendicularFoot(pt0, line0).get_center()
    return Line(perp_foot, perp_unit + perp_foot)

class Test(GraphScene, MovingCameraScene):
    CONFIG = {
        "x_min": -10, "x_max": 10, "x_axis_width": 20,
        "y_min": -5, "y_max": 5, "y_axis_height": 10,
        "axes_color": BLACK,
        "graph_origin": ORIGIN,
    }
    def setup(self):
        GraphScene.setup(self)
        MovingCameraScene.setup(self)
    
    def construct(self):
        self.setup_axes()
        self.play(FadeOut(self.axes))
        self.wait()
        nplane = VGroup(DoubleArrow(8 * LEFT, 8 * RIGHT), DoubleArrow(8 * DOWN, 8 * UP)).set_color(GREY)
        self.play(Write(nplane))
        self.wait()
        planegrid = VGroup()
        for i in range(-10, 10 + 1):
            planegrid.add(Line(8 * UP + i * RIGHT, 8 * DOWN + i * RIGHT))
        for i in range(-5, 5 + 1):
            planegrid.add(Line(i * UP + 10 * RIGHT, i * UP + 10 * LEFT))
        planegrid.set_color(GREY)
        planegrid.fade(0.875)
        self.play(ShowCreation(planegrid), run_time = 3)
        self.wait()
        parab = ParametricFunction(lambda t: np.array([t, np.sqrt(t), 0]), t_min = 0, t_max = 8, color = BLUE)
        self.play(
            Write(parab),
            #self.camera_frame.shift, 1 * UP
        )
        self.play(parab.flip, RIGHT, {"about_point": ORIGIN})
        self.wait()
        semicub = ParametricFunction(lambda t: np.array([t, 2 * t * np.sqrt(t) / 3, 0]), t_min = 0, t_max = 8, color = BLUE)
        semicub.set_color(GREEN)
        self.play(Write(semicub))
        self.play(
            self.camera_frame.scale, 0.75
        )
        def grpcreator(x0):
            x = x0
            c = Circle().shift(x * RIGHT)
            vl = Line(10 * UP, 10 * DOWN).shift(x * RIGHT)
            sl1 = Line(x * RIGHT + RIGHT, x * RIGHT + np.sqrt(x) * DOWN)
            templen = sl1.get_length()
            sl1 = Line(sl1.get_start() + 10 * (sl1.get_end() - sl1.get_start()) / templen, sl1.get_start() - 10 * (sl1.get_end() - sl1.get_start()) / templen)
            sl2 = ParallalLine(Dot().shift(x * RIGHT + 2 * x * np.sqrt(x) * UP / 3), sl1)
            sl2 = Line(sl2.get_start() + 10 * (sl2.get_end() - sl2.get_start()), sl2.get_start() - 10 * (sl2.get_end() - sl2.get_start()))
            sl3 = ParallalLine(Dot(c.get_center()), sl1)
            #sl3 = ParallalLine(Dot().shift(x * RIGHT), sl1)
            sl3.rotate(90 * DEGREES, about_point = c.get_center())
            sl4 = sl2.copy().move_to(sl3.get_end())
            return VGroup(c, vl, sl1, sl2, sl3, sl4)
        self.play(Write(grpcreator(2)))
        self.wait(5)

class Intro(GraphScene, MovingCameraScene):
    CONFIG = {
        "x_min": -10, "x_max": 10, "x_axis_width": 20,
        "y_min": -5, "y_max": 5, "y_axis_height": 10,
        "axes_color": BLACK,
        "graph_origin": ORIGIN,
    }
    def setup(self):
        GraphScene.setup(self)
        MovingCameraScene.setup(self)
    
    def construct(self):
        self.setup_axes()
        self.play(FadeOut(self.axes))
        self.wait()
        nplane = VGroup(DoubleArrow(8 * LEFT, 8 * RIGHT), DoubleArrow(5.5 * DOWN, 5.5 * UP)).set_color(GREY)
        self.play(Write(nplane))
        self.wait()
        planegrid = VGroup()
        for i in range(-10, 10 + 1):
            planegrid.add(Line(5 * UP + i * RIGHT, 5 * DOWN + i * RIGHT))
        for i in range(-5, 5 + 1):
            planegrid.add(Line(i * UP + 10 * RIGHT, i * UP + 10 * LEFT))
        planegrid.set_color(GREY)
        planegrid.fade(0.875)
        self.play(ShowCreation(planegrid), run_time = 3)
        def cyc_generate(rtime = 10, colorpref = GREEN):
            t = ValueTracker(-3 * PI)
            c = Circle()
            l = Line(ORIGIN, DOWN)
            cg = VGroup(c, l)
            cg.rotate_about_origin(-t.get_value())
            cg.shift(t.get_value() * RIGHT + UP)
            path = VMobject(color = colorpref)
            pathd = Dot(cg[-1].get_end(), radius = DEFAULT_SMALL_DOT_RADIUS)
            path.set_points_as_corners([pathd.get_center(), pathd.get_center() + 0.001 * UP])
            pathgp = VGroup(path, pathd)
            def movcircle(obj):
                tc, tl = obj
                tempg = VGroup(Circle(), Line(ORIGIN, DOWN))
                tempg.rotate_about_origin(-t.get_value())
                tempg.shift(t.get_value() * RIGHT + UP)
                obj.become(tempg)
            def pathupdater(obj):
                pat, patd = obj
                patd.move_to(cg[-1].get_end())
                oldpath = pat.copy()
                oldpath.append_vectorized_mobject(Line(oldpath.points[-1], patd.get_center()))
                pat.become(oldpath)
            self.play(Write(cg), Write(pathgp))
            cg.add_updater(movcircle)
            pathgp.add_updater(pathupdater)
            self.add(cg, pathgp)
            self.play(t.set_value, 6 * PI, run_time = rtime, rate_func = linear)
            cg.clear_updaters()
            pathgp.clear_updaters()
            self.play(FadeOut(cg), FadeOut(pathd))
            self.remove(path)
            return path
        '''logspi = ParametricFunction(lambda t: np.array([np.exp(0.2 * t) * np.cos(t), np.exp(0.2 * t) * np.sin(t), 0]), t_min = -4 * PI, t_max = 4 * PI, color = BLUE)
        dot = Dot(color = BLUE)
        dot.move_to(logspi.get_start())
        dot.add_updater(lambda m: m.move_to(logspi.get_end()))
        sline = Line(ORIGIN, dot.get_center(), color = GREEN)
        sline.add_updater(lambda m: m.become(Line(ORIGIN, dot.get_center(), color = GREEN)))
        lline = Line(ORIGIN, 10 * RIGHT).fade(0.75)
        lline.add_updater(lambda m: m.become(Line(ORIGIN, 10 * (sline.get_end() - sline.get_start()) / sline.get_length()).fade(0.75)))
        tracinggrp = VGroup(dot, sline, lline)
        self.play(Write(tracinggrp))
        self.wait()
        self.play(ShowCreation(logspi), run_time = 10, rate_func = linear)'''
        def logspigenerate(rtime = 10, colorpref = BLUE):
            spiang = 80 * DEGREES
            omega = 1
            def grp_creator(temp):
                lline = Line(ORIGIN, 10 * RIGHT).fade(0.75)
                leng = np.exp(omega * temp / np.tan(spiang))
                sline = Line(ORIGIN, leng * RIGHT).set_color(GREEN)
                tgrp = VGroup(lline, sline)
                tgrp.rotate_about_origin(omega * temp)
                return tgrp
            t = ValueTracker(-6 * PI)
            cg = grp_creator(t.get_value())
            path = VMobject(color = colorpref)
            pathd = Dot(cg[-1].get_end(), color = RED)
            path.set_points_as_corners([pathd.get_center(), pathd.get_center() + 0.001 * UP])
            pathgp = VGroup(path, pathd)
            def grp_updater(obj):
                temp = grp_creator(t.get_value())
                obj.become(temp)
            def path_updater(obj):
                pat, patd = obj
                patd.move_to(cg[-1].get_end())
                oldpath = pat.copy()
                oldpath.append_vectorized_mobject(Line(oldpath.points[-1], patd.get_center()))
                oldpath.make_smooth()
                pat.become(oldpath)
            self.play(
                t.set_value, 4 * PI,
                UpdateFromFunc(cg, grp_updater),
                UpdateFromFunc(pathgp, path_updater),
                run_time = rtime,
                rate_func = linear,
            )
            self.play(FadeOut(cg), FadeOut(pathd))
            self.remove(path)
            return path
        logspi = logspigenerate()
        self.add(logspi)
        cycl = cyc_generate()
        self.add(cycl)
        self.wait(5)

class NotionOfLength(GraphScene, MovingCameraScene):
    CONFIG = {
        "x_min": -10, "x_max": 10, "x_axis_width": 20,
        "y_min": -5, "y_max": 5, "y_axis_height": 10,
        "axes_color": BLACK,
        "graph_origin": ORIGIN,
    }
    def setup(self):
        GraphScene.setup(self)
        MovingCameraScene.setup(self)
    
    def construct(self):
        self.setup_axes()
        self.play(FadeOut(self.axes))
        self.wait()
        stline = Line(5 * LEFT + 3 * DOWN, 2 * UP + 2 * LEFT, color = YELLOW)
        curv = VMobject(color = YELLOW)
        curv.set_points_smoothly([1 * RIGHT + 2 * UP, 4 * RIGHT + 0 * UP, 0.5 * RIGHT + 1 * DOWN, 3 * DOWN + 3 * RIGHT])
        curv.shift(RIGHT)
        llabel = TexMobject("\\underline{\\text{Straight line}}")
        clabel = TexMobject("\\underline{\\text{'squiggly' Curve}}")
        llabel.shift(3 * UP + 3 * LEFT)
        clabel.shift(3 * UP + 3 * RIGHT)
        self.play(Write(llabel), Write(stline))
        self.wait()
        self.play(
            ReplacementTransform(stline.copy(), curv),
            Write(clabel))
        self.wait()
        lenlabel = TexMobject("\\Leftarrow\\text{ Length }\\Rightarrow")
        self.play(Write(lenlabel))
        self.wait()
        self.play(*[FadeOut(obj) for obj in [llabel, clabel, stline, curv, lenlabel]])
        self.wait()
        title = TexMobject("\\underline{\\text{Descartes' assertion}}").shift(3.5 * UP)
        self.play(Write(title))
        self.wait()
        mechgrp = VGroup(
            TextMobject("Mechanical Curves"),
            TextMobject("defined by the motion of simple objects").scale(0.875),
            TextMobject("e.g.: Logarithmic Spiral, Cycloid").scale(0.875)
        )
        mechgrp[0].set_color(YELLOW)
        mechgrp.arrange(direction = DOWN, buff = 0.375)
        for obj in mechgrp[1:]:
            obj.align_to(mechgrp[0], LEFT)
        for obj in mechgrp[1:]:
            obj.shift(0.5 * RIGHT)
        mechgrp.align_on_border(LEFT)
        mechgrp.shift(UP + 2 * RIGHT)
        self.play(Write(mechgrp))
        self.wait()
        mathgrp = VGroup(
            TextMobject("Mathematical Curves"),
            TextMobject("defined by mathematical equations").scale(0.875),
            TextMobject("e.g.: Power functions, Logarithms").scale(0.875)
        )
        mathgrp[0].set_color(YELLOW)
        mathgrp.arrange(direction = DOWN, buff = 0.375)
        for obj in mathgrp[1:]:
            obj.align_to(mathgrp[0], LEFT)
        for obj in mathgrp[1:]:
            obj.shift(0.5 * RIGHT)
        mathgrp.align_on_border(LEFT)
        mathgrp.shift(2 * DOWN + 2 * RIGHT)
        self.play(Write(mathgrp))
        self.wait()
        cr = Cross(mechgrp)
        cr.rotate_in_place(-3 * DEGREES)
        self.play(Write(cr))
        self.wait()
        self.play(
            *[FadeOut(obj) for obj in [cr, mechgrp, mathgrp, title]]
        )
        self.wait()
        neilname = VGroup(
            TextMobject("William Neil"),
            TextMobject("(1637-1670)")
        )
        neilname.arrange(DOWN)
        self.play(Write(neilname))
        self.wait(5)

class Semicubical(GraphScene, MovingCameraScene):
    CONFIG = {
        "x_min": -10, "x_max": 10, "x_axis_width": 20,
        "y_min": -5, "y_max": 5, "y_axis_height": 10,
        "axes_color": BLACK,
        "graph_origin": ORIGIN,
    }
    def setup(self):
        GraphScene.setup(self)
        MovingCameraScene.setup(self)
    
    def construct(self):
        self.setup_axes()
        self.play(FadeOut(self.axes))
        self.wait()
        nplane = VGroup(DoubleArrow(8 * LEFT, 8 * RIGHT), DoubleArrow(5.5 * DOWN, 5.5 * UP)).set_color(GREY)
        self.play(Write(nplane))
        self.wait()
        planegrid = VGroup()
        for i in range(-10, 10 + 1):
            planegrid.add(Line(5 * UP + i * RIGHT, 5 * DOWN + i * RIGHT))
        for i in range(-5, 5 + 1):
            planegrid.add(Line(i * UP + 10 * RIGHT, i * UP + 10 * LEFT))
        planegrid.set_color(GREY)
        planegrid.fade(0.875)
        self.play(ShowCreation(planegrid), run_time = 3)
        self.wait()
        scub = ParametricFunction(lambda t: np.array([t, t * np.sqrt(t), 0]), t_min = 0, t_max = 5, color = BLUE)
        scub1 = ParametricFunction(lambda t: np.array([t, -t * np.sqrt(t), 0]), t_min = 0, t_max = 5, color = BLUE)
        neilname = VGroup(
            TextMobject("William Neil"),
            TextMobject("(1637-1670)")
        )
        neilname.arrange(DOWN)
        neilname.add_background_rectangle()
        self.play(Write(neilname))
        self.wait()
        self.play(
            neilname.shift, 4 * LEFT + 2 * UP,
            ShowCreation(scub),
            ShowCreation(scub1),
        )
        self.wait()
        semieqn = TexMobject("y=\\pm a x^{3/2}")
        semieqn.add_background_rectangle()
        semieqn.move_to(3.5 * RIGHT + UP + semieqn.get_center() - semieqn.get_critical_point(DOWN))
        self.play(
            Write(semieqn),
            FadeOut(neilname))
        self.wait(5)

class NeilProof(GraphScene, MovingCameraScene):
    CONFIG = {
        "x_min": -10, "x_max": 10, "x_axis_width": 20,
        "y_min": -5, "y_max": 5, "y_axis_height": 10,
        "axes_color": BLACK,
        "graph_origin": ORIGIN,
    }
    def setup(self):
        GraphScene.setup(self)
        MovingCameraScene.setup(self)
    
    def construct(self):
        self.setup_axes()
        self.play(FadeOut(self.axes))
        self.wait()
        nplane = VGroup(DoubleArrow(8 * LEFT, 10 * RIGHT), DoubleArrow(5.5 * DOWN, 5.5 * UP)).set_color(GREY)
        self.play(Write(nplane))
        self.wait()
        planegrid = VGroup()
        for i in range(-10, 10 + 1):
            planegrid.add(Line(5 * UP + i * RIGHT, 5 * DOWN + i * RIGHT))
        for i in range(-5, 5 + 1):
            planegrid.add(Line(i * UP + 10 * RIGHT, i * UP + 10 * LEFT))
        planegrid.set_color(GREY)
        planegrid.fade(0.875)
        self.play(ShowCreation(planegrid), run_time = 3)
        self.wait()
        parab = self.get_graph(
            lambda x: x ** 2,
            x_min = -3,
            x_max = 3,
            color = BLUE,
        )
        self.play(Write(parab))
        self.wait()
        self.play(
            parab.flip, UR, {"about_point": ORIGIN},
            self.camera_frame.shift, UP + RIGHT)
        self.wait()
        parabeqn = TexMobject("y^2", "=", "x")
        parabeqn.add_background_rectangle()
        parabeqn.move_to(4 * RIGHT + DOWN)
        self.play(Write(parabeqn))
        self.wait()
        hparab = self.get_graph(
            lambda x: x ** 2,
            x_min = 0,
            x_max = 3,
            color = BLUE,
        ).flip(UR, about_point = ORIGIN)
        hparabeqn = TexMobject("y", "=", "\\sqrt{x}")
        hparabeqn.add_background_rectangle()
        hparabeqn.move_to(parabeqn[1].get_center() + hparabeqn.get_center() - hparabeqn[1].get_center())
        self.add(hparab)
        self.play(
            FadeOut(parab),
            Transform(parabeqn, hparabeqn)
        )
        self.wait()
        def area_creator(temp):
            narea = VMobject(fill_color = GREEN, fill_opacity = 0.75, stroke_width = 1)
            nareapts = [np.array([i, np.sqrt(i), 0]) for i in np.arange(0, temp + 0.1, 0.1)]
            #nareapts = [ORIGIN] + nareapts + [temp * RIGHT + 2 * temp * np.sqrt(temp) * UP / 3, temp * RIGHT]
            nareapts[-1][0] = temp
            nareapts[-1][1] = np.sqrt(temp)
            nareapts = [ORIGIN] + nareapts + [temp * RIGHT]
            #print(nareapts[-1][0], nareapts[-2][0])
            #nareapts.append(temp * RIGHT)
            narea.set_points_as_corners(nareapts)
            #nrectan = Rectangle(width = temp, height = np.sqrt(temp)).shift(temp * RIGHT / 2 + np.sqrt(temp) * UP / 2).set_color(PINK)
            return narea
        xval = 3
        area = area_creator(xval)
        self.play(Write(area))
        self.wait()
        rectan = Rectangle(width = xval, height = np.sqrt(xval)).shift(xval * RIGHT / 2 + np.sqrt(xval) * UP / 2).set_color(GREEN)
        self.play(ShowCreation(rectan))
        self.wait()
        braces = VGroup(
            Brace(Line(ORIGIN, xval * RIGHT), DOWN),
            Brace(Line(xval * RIGHT, xval * RIGHT + np.sqrt(xval) * UP), RIGHT)
        )
        braces.add(TexMobject("x").next_to(braces[0], DOWN))
        braces.add(TexMobject("\\sqrt{x}").next_to(braces[1], RIGHT))
        self.play(
            Write(braces),
            parabeqn.move_to, 6.5 * RIGHT + 0.5 * UP)
        self.wait()
        arealabel = TexMobject("A(x)", "=", "\\frac{2}{3}x^{3/2}")
        arealabel.move_to(3 * LEFT + UP + arealabel.get_center() - arealabel[1].get_center())
        self.play(Write(arealabel))
        self.wait()
        self.play(
            *[FadeOut(obj) for obj in [rectan, area]],
            *[FadeOutAndShiftDown(obj) for obj in [braces[0], braces[2]]],
            *[FadeOutAndShift(obj, RIGHT) for obj in [parabeqn, braces[1], braces[3]]],
        )
        self.wait()
        t = ValueTracker(0.1)
        def pgrp_creator(temp):
            ar = area_creator(temp)
            l = Line(ORIGIN, 2 * temp * np.sqrt(temp) * UP / 3).fade(0.5)
            l.shift(temp * RIGHT)
            return VGroup(ar, l)
        def pgrpupdater(obj):
            temp = pgrp_creator(t.get_value())
            obj.become(temp)
        argrp = pgrp_creator(t.get_value())
        path = VMobject(color = BLUE)
        pathd = Dot(argrp[-1].get_end(), radius = DEFAULT_SMALL_DOT_RADIUS)
        path.set_points_as_corners([pathd.get_center(), pathd.get_center() + 0.001 * UP])
        pathgp = VGroup(path, pathd)
        def pathupdater(obj):
            pat, patd = obj
            patd.move_to(argrp[-1].get_end())
            oldpath = pat.copy()
            oldpath.append_vectorized_mobject(Line(oldpath.points[-1], patd.get_center()))
            pat.become(oldpath)
        self.play(
            t.set_value, 5,
            UpdateFromFunc(argrp, pgrpupdater),
            UpdateFromFunc(pathgp, pathupdater),
            run_time = 6
        )
        self.wait()
        spara = ParametricFunction(lambda t: np.array([t, 2 * t * np.sqrt(t) / 3, 0]), t_min = 0, t_max = 5, color = BLUE)
        self.add(spara)
        self.play(
            FadeOut(pathgp),
            FadeOut(argrp)
        )
        self.wait()
        self.play(
            hparab.flip, RIGHT, {"about_point": ORIGIN},
            spara.set_color, PURPLE,
        )
        self.wait(5)

class TangentIntroduction(GraphScene, MovingCameraScene):
    CONFIG = {
        "x_min": -10, "x_max": 10, "x_axis_width": 20,
        "y_min": -5, "y_max": 5, "y_axis_height": 10,
        "axes_color": BLACK,
        "graph_origin": ORIGIN,
    }
    def setup(self):
        GraphScene.setup(self)
        MovingCameraScene.setup(self)
    
    def construct(self):
        self.setup_axes()
        self.play(FadeOut(self.axes))
        self.wait()
        nplane = VGroup(DoubleArrow(8 * LEFT, 10 * RIGHT), DoubleArrow(5.5 * DOWN, 5.5 * UP)).set_color(GREY)
        #self.play(Write(nplane))
        #self.wait()
        planegrid = VGroup()
        for i in range(-10, 10 + 1):
            planegrid.add(Line(5 * UP + i * RIGHT, 5 * DOWN + i * RIGHT))
        for i in range(-5, 5 + 1):
            planegrid.add(Line(i * UP + 10 * RIGHT, i * UP + 10 * LEFT))
        planegrid.set_color(GREY)
        planegrid.fade(0.875)
        #self.add(planegrid)
        #self.play(ShowCreation(planegrid), run_time = 3)
        spara = ParametricFunction(lambda t: np.array([t, 2 * t * np.sqrt(t) / 3, 0]), t_min = 0, t_max = 5, color = PURPLE)
        arealabel = TexMobject("A(x)", "=", "\\frac{2}{3}x^{3/2}")
        arealabel.move_to(3 * LEFT + UP + arealabel.get_center() - arealabel[1].get_center())
        hparab = self.get_graph(
            lambda x: x ** 2,
            x_min = 0,
            x_max = 3,
            color = BLUE,
        ).flip(UR, about_point = ORIGIN)
        hparab.flip(RIGHT, about_point = ORIGIN)
        self.camera_frame.shift(UP + RIGHT)
        self.add(nplane, planegrid, hparab, spara, arealabel)
        self.wait()
        t = ValueTracker(0)
        tangentpoint = Dot(color = RED)
        tangentline = Line(12 * LEFT, 12 * RIGHT, color = RED).save_state()
        self.play(Write(tangentpoint), Write(tangentline))
        tangentpoint.add_updater(lambda m: m.move_to(t.get_value() * RIGHT + 2 * t.get_value() * np.sqrt(t.get_value()) * UP / 3))
        tangentline.add_updater(lambda m: m.restore().rotate_about_origin(np.arctan(np.sqrt(t.get_value()))).move_to(t.get_value() * RIGHT + 2 * t.get_value() * np.sqrt(t.get_value()) * UP / 3))
        self.play(t.set_value, 1, run_time = 1.5)
        self.play(t.set_value, 2, run_time = 1.5)
        self.play(t.set_value, 3, run_time = 1.5)
        self.wait(0.5)
        self.play(t.set_value, 2, run_time = 1.5)
        self.wait()
        pascal = ImageMobject("Blaise_Pascal_Versailles").scale(2.5).shift(UR).shift(3 * LEFT)
        fermat = ImageMobject("Pierre_de_Fermat").scale(2.5).shift(UR).shift(3 * RIGHT)
        barrow = ImageMobject("Isaac_Barrow_by_Mary_Beale").scale(3).shift(UR).shift(0.5 * UP)
        pascalname = TextMobject("Blaise Pascal").add_background_rectangle().next_to(pascal, DOWN)
        fermatname = TextMobject("Pierre de Fermat").add_background_rectangle().next_to(fermat, DOWN)
        barrowname = TextMobject("Isaac Barrow").add_background_rectangle().next_to(barrow, DOWN)
        self.play(
            *[FadeInFrom(obj, LEFT) for obj in [pascal, pascalname]],
            *[FadeInFrom(obj, RIGHT) for obj in [fermat, fermatname]],
        )
        self.wait()
        self.play(
            *[FadeOutAndShift(obj, LEFT) for obj in [pascal, pascalname]],
            *[FadeOutAndShift(obj, RIGHT) for obj in [fermat, fermatname]],
            FadeInFromLarge(barrow),
            Write(barrowname)
        )
        self.wait()
        self.play(FadeOut(barrow), FadeOut(barrowname))
        self.play(FadeOut(tangentpoint), FadeOut(tangentline))
        self.wait(5)

class TangentConstruction(GraphScene, MovingCameraScene):
    CONFIG = {
        "x_min": -10, "x_max": 10, "x_axis_width": 20,
        "y_min": -5, "y_max": 5, "y_axis_height": 10,
        "axes_color": BLACK,
        "graph_origin": ORIGIN,
    }
    def setup(self):
        GraphScene.setup(self)
        MovingCameraScene.setup(self)
    
    def construct(self):
        self.setup_axes()
        self.play(FadeOut(self.axes))
        self.wait()
        nplane = VGroup(DoubleArrow(8 * LEFT, 10 * RIGHT), DoubleArrow(5.5 * DOWN, 5.5 * UP)).set_color(GREY)
        #self.play(Write(nplane))
        #self.wait()
        planegrid = VGroup()
        for i in range(-10, 10 + 1):
            planegrid.add(Line(5 * UP + i * RIGHT, 5 * DOWN + i * RIGHT))
        for i in range(-5, 5 + 1):
            planegrid.add(Line(i * UP + 10 * RIGHT, i * UP + 10 * LEFT))
        planegrid.set_color(GREY)
        planegrid.fade(0.875)
        #self.add(planegrid)
        #self.play(ShowCreation(planegrid), run_time = 3)
        spara = ParametricFunction(lambda t: np.array([t, 2 * t * np.sqrt(t) / 3, 0]), t_min = 0, t_max = 5, color = PURPLE)
        arealabel = TexMobject("A(x)", "=", "\\frac{2}{3}x^{3/2}")
        arealabel.move_to(3 * LEFT + UP + arealabel.get_center() - arealabel[1].get_center())
        hparab = self.get_graph(
            lambda x: x ** 2,
            x_min = 0,
            x_max = 3,
            color = BLUE,
        ).flip(UR, about_point = ORIGIN)
        hparab.flip(RIGHT, about_point = ORIGIN)
        self.camera_frame.shift(UP + RIGHT)
        self.add(nplane, planegrid, hparab, spara, arealabel)
        t = ValueTracker(2.5)
        xdot = Dot()
        xdot.add_updater(lambda m: m.move_to(t.get_value() * RIGHT))
        xlab = TexMobject("O").scale(0.75).add_background_rectangle()
        xlab.add_updater(lambda m: m.next_to(xdot, DR, buff = 0.0))
        line = Line(xdot.get_center(), xdot.get_center() + np.power(t.get_value(), 1.5) * 2 * UP / 3).fade(0.5)
        line.add_updater(lambda m: m.become(Line(xdot.get_center(), xdot.get_center() + np.power(t.get_value(), 1.5) * 2 * UP / 3).fade(0.5)))
        pdot = Dot(line.get_end())
        pdot.add_updater(lambda m: m.move_to(line.get_end()))
        plab = TexMobject("P").scale(0.75).add_background_rectangle()
        plab.add_updater(lambda m: m.next_to(pdot, UL, buff = 0.0))
        self.play(
            *[Write(obj) for obj in [xdot, xlab, line, pdot, plab]]
        )
        vline = Line(10 * UP, 10 * DOWN, color = GREEN).save_state()
        vline.add_updater(lambda m: m.restore().shift(t.get_value() * RIGHT))
        pd = Dot()
        pd.add_updater(lambda m: m.move_to(t.get_value() * RIGHT + np.sqrt(t.get_value()) * DOWN))
        pddot = TexMobject("P'").scale(0.75).add_background_rectangle()
        pddot.add_updater(lambda m: m.next_to(pd, DR, buff = 0))
        self.play(Write(vline), Write(pd), Write(pddot))
        self.wait()
        c = Circle()
        c.add_updater(lambda m: m.move_to(xdot.get_center()))
        qdot = xdot.copy().shift(RIGHT)
        qdot.add_updater(lambda m: m.move_to(xdot.get_center() + RIGHT))
        #qddot = xdot.copy().shift(LEFT)
        #qddot.add_updater(lambda m: m.move_to(xdot.get_center() + LEFT))
        qlab = TexMobject("Q").scale(0.75).add_background_rectangle()
        #qdlab = TexMobject("Q'").scale(0.75).add_background_rectangle()
        qlab.add_updater(lambda m: m.next_to(qdot, DOWN, buff = 0.05))
        #qdlab.add_updater(lambda m: m.next_to(qddot, DOWN, buff = 0.05))
        self.play(
            #*[Write(obj) for obj in [c, qdot, qddot, qlab, qdlab]]
            *[Write(obj) for obj in [c, qdot, qlab]]
        )
        self.wait()
        def lineextend(s0, e0, color = RED):
            s, e = s0, e0
            temp = Line(s0, e0)
            leng = temp.get_length()
            return Line(s0 + 100 * (e0 - s0) / leng, s0 - 100 * (e0 - s0) / leng, color = color)
        qpdline = Line()
        qpdline.add_updater(lambda m: m.become(lineextend(qdot.get_center(), pd.get_center(), color = ORANGE)))
        self.play(Write(qpdline))
        self.wait()
        tangentline = lineextend(qdot.get_center(), pd.get_center(), color = ORANGE)
        self.play(tangentline.move_to, pdot.get_center())
        tangentline.add_updater(lambda m: m.become(lineextend(qdot.get_center(), pd.get_center(), color = ORANGE).move_to(pdot.get_center())))
        self.wait()
        self.play(t.increment_value, 1, run_time = 1.5)
        self.wait(0.5)
        self.play(t.increment_value, -2, run_time = 3)
        self.wait(0.5)
        self.play(t.increment_value, 1, run_time = 1.5)
        self.wait()
        def ctangcreator(qob, pdob, xob):
            templine = Line(pdob, qob)
            templine = Line(templine.get_start(), templine.get_start() + (templine.get_end() - templine.get_start()) / templine.get_length())
            templine.move_to(xob + templine.get_center() - templine.get_start())
            templine.rotate(90 * DEGREES, about_point = xob)
            templine.rotate(90 * DEGREES, about_point = templine.get_end())
            templine = lineextend(templine.get_start(), templine.get_end(), color = ORANGE)
            return templine
        ctang = ctangcreator(qdot.get_center(), pd.get_center(), xdot.get_center())
        qpdline.clear_updaters()
        self.play(Write(ctang), FadeOut(qpdline))
        self.wait()
        ctang.add_updater(lambda m: m.become(ctangcreator(qdot.get_center(), pd.get_center(), xdot.get_center())))
        interdot = Dot()
        interdot.add_updater(lambda m: m.move_to(line_intersection(ctang.get_start_and_end(), vline.get_start_and_end())))
        tlab = TexMobject("T").scale(0.75).add_background_rectangle()
        tlab.add_updater(lambda m: m.next_to(interdot, DR, buff = 0))
        self.play(Write(interdot), Write(tlab))
        self.wait()
        self.play(t.increment_value, 1, run_time = 1.5)
        self.wait(0.5)
        self.play(
            arealabel.shift, UP,
            t.set_value, 1e-5, run_time = 3)
        path = VMobject(color = YELLOW)
        pathd = Dot(interdot.get_center(), radius = DEFAULT_SMALL_DOT_RADIUS)
        path.set_points_as_corners([pathd.get_center(), pathd.get_center() + 0.001 * UP])
        pathgp = VGroup(path, pathd)
        def pathupdater(obj):
            pat, patd = obj
            patd.move_to(interdot.get_center())
            oldpath = pat.copy()
            oldpath.append_vectorized_mobject(Line(oldpath.points[-1], patd.get_center()))
            pat.become(oldpath)
        self.play(Write(pathgp))
        pathgp.add_updater(pathupdater)
        self.add(pathgp)
        self.play(t.set_value, 10, run_time = 7.5)
        self.wait()
        pathgp.clear_updaters()
        auxcur = ParametricFunction(lambda t: np.array([t, np.sqrt(1 + t), 0]), t_min = 0, t_max = 10, color = YELLOW)
        self.add(auxcur)
        self.remove(path)
        self.play(t.set_value, 3, run_time = 5)
        self.wait(5)

class AuxillaryCurve(GraphScene, MovingCameraScene):
    CONFIG = {
        "x_min": -10, "x_max": 10, "x_axis_width": 20,
        "y_min": -5, "y_max": 5, "y_axis_height": 10,
        "axes_color": BLACK,
        "graph_origin": ORIGIN,
    }
    def setup(self):
        GraphScene.setup(self)
        MovingCameraScene.setup(self)
    
    def construct(self):
        self.setup_axes()
        self.play(FadeOut(self.axes))
        self.wait()
        nplane = VGroup(DoubleArrow(8 * LEFT, 10 * RIGHT), DoubleArrow(5.5 * DOWN, 5.5 * UP)).set_color(GREY)
        planegrid = VGroup()
        for i in range(-10, 10 + 1):
            planegrid.add(Line(5 * UP + i * RIGHT, 5 * DOWN + i * RIGHT))
        for i in range(-5, 5 + 1):
            planegrid.add(Line(i * UP + 10 * RIGHT, i * UP + 10 * LEFT))
        planegrid.set_color(GREY)
        planegrid.fade(0.875)
        spara = ParametricFunction(lambda t: np.array([t, 2 * t * np.sqrt(t) / 3, 0]), t_min = 0, t_max = 5, color = PURPLE)
        arealabel = TexMobject("A(x)", "=", "\\frac{2}{3}x^{3/2}")
        arealabel.move_to(3 * LEFT + UP + arealabel.get_center() - arealabel[1].get_center())
        arealabel.shift(UP)
        hparab = self.get_graph(
            lambda x: x ** 2,
            x_min = 0,
            x_max = 3,
            color = BLUE,
        ).flip(UR, about_point = ORIGIN)
        hparab.flip(RIGHT, about_point = ORIGIN)
        self.camera_frame.shift(UP + RIGHT)
        t = ValueTracker(3)
        xdot = Dot(t.get_value() * RIGHT)
        #xdot.add_updater(lambda m: m.move_to(t.get_value() * RIGHT))
        xlab = TexMobject("O").scale(0.75).add_background_rectangle().next_to(xdot, DR, buff = 0.0)
        #xlab.add_updater(lambda m: m.next_to(xdot, DR, buff = 0.0))
        line = Line(xdot.get_center(), xdot.get_center() + np.power(t.get_value(), 1.5) * 2 * UP / 3).fade(0.5)
        #line.add_updater(lambda m: m.become(Line(xdot.get_center(), xdot.get_center() + np.power(t.get_value(), 1.5) * 2 * UP / 3).fade(0.5)))
        pdot = Dot(line.get_end())
        #pdot.add_updater(lambda m: m.move_to(line.get_end()))
        plab = TexMobject("P").scale(0.75).add_background_rectangle().next_to(pdot, UL, buff = 0.0)
        #plab.add_updater(lambda m: m.next_to(pdot, UL, buff = 0.0))
        vline = Line(10 * UP, 10 * DOWN, color = GREEN).shift(t.get_value() * RIGHT)
        pd = Dot(t.get_value() * RIGHT + np.sqrt(t.get_value()) * DOWN)
        pddot = TexMobject("P'").scale(0.75).add_background_rectangle().next_to(pd, DR, buff = 0)
        #pddot.add_updater(lambda m: m.next_to(pd, DR, buff = 0))
        auxcur = ParametricFunction(lambda t: np.array([t, np.sqrt(1 + t), 0]), t_min = 0, t_max = 10, color = YELLOW)
        qdot = xdot.copy().shift(RIGHT)
        #qdot.add_updater(lambda m: m.move_to(xdot.get_center() + RIGHT))
        qlab = TexMobject("Q").scale(0.75).add_background_rectangle().next_to(qdot, DR, buff = 0)
        c = Circle(arc_center = xdot.get_center())
        tdot = Dot(t.get_value() * RIGHT + np.sqrt(1 + t.get_value()) * UP)
        tlab = TexMobject("T").scale(0.75).add_background_rectangle().next_to(tdot, DR, buff = 0)
        tangentline = Line(10 * LEFT, 10 * RIGHT, color = ORANGE).rotate_about_origin(np.arctan(np.sqrt(t.get_value()))).move_to(pdot.get_center())
        ctangentline = Line(10 * LEFT, 10 * RIGHT, color = ORANGE).rotate_about_origin(np.arctan(np.sqrt(t.get_value()))).move_to(tdot.get_center())
        self.add(nplane, planegrid, hparab, spara, arealabel, xdot, xlab, pdot, plab, vline, pd, pddot, auxcur, qdot, qlab, c, tdot, tlab, tangentline, ctangentline)
        self.wait()
        self.play(
            FadeOut(ctangentline),
            FadeOut(auxcur)
        )
        self.wait()
        wid = 0.75
        rectan = Rectangle(length = t.get_value(), width = wid, color = TEAL)
        rectan.move_to(tdot.get_center() + rectan.get_center() - rectan.get_critical_point(UP))
        self.play(ShowCreation(rectan))
        self.wait()
        midline = Line(xdot.get_center(), tdot.get_center())
        lline = midline.copy().shift(wid * LEFT / 2)
        rline = midline.copy().shift(wid * RIGHT / 2)
        linter = line_intersection(lline.get_start_and_end(), tangentline.get_start_and_end())
        rinter = line_intersection(rline.get_start_and_end(), tangentline.get_start_and_end())
        lfad = Line(lline.get_start(), linter).fade(0.75)
        rfad = Line(rline.get_start(), rinter).fade(0.75)
        self.play(ShowCreation(lfad), ShowCreation(rfad))
        self.wait()
        tansl = np.arctan(tangentline.get_slope())
        righttri = VGroup(
            Line(linter, rinter),
            Line(rinter, linter + wid * RIGHT),
            Line(linter + wid * RIGHT, linter),
        ).set_color(TEAL)
        self.play(Write(righttri))
        self.wait()
        fade_val = 7 / 8
        reverse_fade = 1 - 1 / (1 - fade_val)
        radline = PerpendicularLine(xdot, ctangentline).set_color(GOLD)
        tddot = Dot(radline.get_start())
        tdlab = TexMobject("T'").scale(0.75).add_background_rectangle().next_to(tddot.get_center(), UL, buff = 0.05)
        horline = Line(ORIGIN, RIGHT)
        sdot = Dot(line_intersection(horline.get_start_and_end(), ctangentline.get_start_and_end()))
        rdot = Dot(line_intersection(horline.get_start_and_end(), tangentline.get_start_and_end()))
        slab = TexMobject("S").scale(0.75).add_background_rectangle().next_to(sdot, DL, buff = 0.05)
        rlab = TexMobject("R").scale(0.75).add_background_rectangle().next_to(rdot, DL, buff = 0.05)
        rangle = Angle(tdot.get_center(), tddot.get_center(), xdot.get_center(), color = YELLOW, radius = 1 / 4)
        llab = TexMobject("L").scale(0.75).add_background_rectangle().next_to(rinter, RIGHT, buff = 0.05)
        mlab = TexMobject("M").scale(0.75).add_background_rectangle().next_to(linter, LEFT, buff = 0.05)
        nlab = TexMobject("N").scale(0.75).add_background_rectangle().next_to(linter + wid * RIGHT, RIGHT, buff = 0.05)
        self.play(
            rectan.fade, fade_val,
            lfad.fade, fade_val,
            rfad.fade, fade_val,
            #*[FadeIn(obj) for obj in [ctangentline, radline, slab, rlab, tdlab, rdot, sdot, tddot, rangle, llab, mlab, nlab]],
            FadeIn(ctangentline),
        )
        self.wait()
        self.play(
            arealabel.shift, 3.5 * DOWN
        )
        simeqn = TexMobject("\\frac{TO}{OT'}", "=", "\\frac{TS}{SO}", "=", "\\frac{PR}{RO}", "=", "\\frac{LM}{MN}").scale(0.75)
        simeqn.move_to(2 * UL + LEFT)
        self.play(*[FadeIn(obj) for obj in [slab, tddot, tdlab, rangle, radline, sdot]])
        self.play(Write(simeqn[:3]))
        self.wait()
        self.play(
            FadeOut(rangle),
            #*[FadeOut(obj) for obj in [tddot, tdlab, rangle, radline]],
            *[FadeIn(obj) for obj in [rlab, rdot]]
        )
        self.play(Write(simeqn[3:5]))
        self.wait()
        self.play(
            *[FadeOut(obj) for obj in [sdot, slab]],
            *[FadeIn(obj) for obj in [llab, mlab, nlab]]
        )
        self.play(Write(simeqn[5:7]))
        self.wait()
        self.play(
            FadeOut(simeqn[1:3]),
            FadeOut(simeqn[4:6]),
            simeqn[0].next_to, simeqn[3], {"direction": LEFT},
            simeqn[-1].next_to, simeqn[3], {"direction": RIGHT},
            *[FadeOut(obj) for obj in [ctangentline, rlab, rdot]],
            rectan.fade, reverse_fade,
            lfad.fade, reverse_fade,
            rfad.fade, reverse_fade,
        )
        self.wait()
        reeqn = TexMobject("TO", "\\cdot", "MN", "=", "OT'", "\\cdot", "LM").scale(0.75).move_to(2 * UL + LEFT)
        onelab = TexMobject("1").scale(0.75).move_to(reeqn[4].get_center())
        self.play(
            FadeOutAndShiftDown(simeqn[0]),
            FadeOutAndShiftDown(simeqn[-1]),
            FadeOutAndShiftDown(simeqn[3]),
            FadeInFrom(reeqn, UP)
        )
        self.wait()
        self.play(
            *[FadeOut(obj) for obj in [tddot, tdlab, radline]],
            Transform(reeqn[4], onelab)
        )
        self.play(
            FadeOut(reeqn[4]),
            FadeOut(reeqn[5]),
            FadeIn(auxcur),
            reeqn[-1].next_to, reeqn[3], {"direction": RIGHT}
        )
        self.wait()
        self.play(
            *[FadeOut(obj) for obj in [rectan, xdot, xlab, pdot, plab, pd, pddot, c, qdot, qlab, tdot, tlab, llab, mlab, nlab, tangentline, righttri, lfad, rfad, vline]],
            *[FadeOut(obj) for obj in [arealabel, reeqn[:4], reeqn[-1]]]
        )
        self.wait(5)

class AuxRiemann(GraphScene, MovingCameraScene):
    CONFIG = {
        "x_min": -10, "x_max": 10, "x_axis_width": 20,
        "y_min": -5, "y_max": 5, "y_axis_height": 10,
        "axes_color": BLACK,
        "graph_origin": ORIGIN,
    }
    def setup(self):
        GraphScene.setup(self)
        MovingCameraScene.setup(self)
    
    def construct(self):
        self.setup_axes()
        self.play(FadeOut(self.axes))
        self.wait()
        nplane = VGroup(DoubleArrow(8 * LEFT, 15 * RIGHT), DoubleArrow(5.5 * DOWN, 8 * UP)).set_color(GREY)
        planegrid = VGroup()
        for i in range(-10, 15 + 1):
            planegrid.add(Line(8 * UP + i * RIGHT, 5 * DOWN + i * RIGHT))
        for i in range(-5, 8 + 1):
            planegrid.add(Line(i * UP + 15 * RIGHT, i * UP + 10 * LEFT))
        planegrid.set_color(GREY)
        planegrid.fade(0.875)
        spara = ParametricFunction(lambda t: np.array([t, 2 * t * np.sqrt(t) / 3, 0]), t_min = 0, t_max = 5, color = PURPLE)
        arealabel = TexMobject("A(x)", "=", "\\frac{2}{3}x^{3/2}")
        arealabel.move_to(3 * LEFT + UP + arealabel.get_center() - arealabel[1].get_center())
        arealabel.shift(UP)
        hparab = self.get_graph(
            lambda x: x ** 2,
            x_min = 0,
            x_max = 3,
            color = BLUE,
        ).flip(UR, about_point = ORIGIN)
        hparab.flip(RIGHT, about_point = ORIGIN)
        self.camera_frame.shift(UP + RIGHT)
        scub = self.get_graph(
            lambda x: 2 * np.sqrt(x ** 3) / 3,
            x_min = 0,
            x_max = 10,
            color = PURPLE,
        )
        auxcur = self.get_graph(
            lambda x: np.sqrt(1 + x),
            x_min = 0,
            x_max = 15,
            color = YELLOW,
        )
        self.add(nplane, planegrid, hparab, scub, auxcur)
        def rectanglesandsegments(minx = 0.0, maxx = 10.0, dx = 0.1, fill_opacity = 0.95):
            rectangles = VGroup()
            segments = VGroup()
            x_range = np.arange(minx, maxx, dx)
            colors = color_gradient([BLUE, GREEN], len(x_range))
            for x, color in zip(x_range, colors):
                sample_point = x + 0.5 * dx
                y_sample = 2 * np.sqrt(sample_point ** 3) / 3
                tansl = np.sqrt(sample_point)
                graph_point = self.input_to_graph_point(sample_point, auxcur)
                seg = Line(x * RIGHT, (x + dx) * RIGHT + dx * tansl * UP).move_to(sample_point * RIGHT + y_sample * UP).set_color(GREEN)
                points = VGroup(*list(map(VectorizedPoint, [
                    self.coords_to_point(x, 0),
                    self.coords_to_point(x + 1.001 * dx, 0),
                    graph_point
                ])))

                rect = Rectangle()
                rect.replace(points, stretch=True)
                #rect.set_fill(color, opacity = fill_opacity)
                rect.set_stroke(TEAL)
                rectangles.add(rect)
                segments.add(seg)
            return rectangles, segments
        def getthelist(n_iter = 6, max_dx = 2):
            rectlist, seglist = [], []
            for n in range(n_iter):
                tempr, temps = rectanglesandsegments(dx = float(max_dx) / 2 ** n)
                rectlist += [tempr]
                seglist += [temps]
            return rectlist, seglist
        rrect, rsegs = getthelist()
        nrrect, nrsegs = rrect[0], rsegs[0]
        self.play(
            self.camera_frame.shift, 2 * RIGHT + UP,
            #self.camera_frame.scale, 1.25
        )
        self.wait()
        self.play(Write(nrrect), Write(nrsegs))
        self.wait()
        for i in range(1, len(rrect)):
            self.play(
                Transform(nrrect, rrect[i], lag_ratio = 0.25),
                Transform(nrsegs, rsegs[i], lag_ratio = 0.25)
            )
            self.wait(0.5)
        self.wait()
        self.play(
            FadeOut(nrrect),
            FadeOut(nrsegs),
            self.camera_frame.shift, 2 * LEFT + DOWN,
        )
        self.wait(5)

class AuxEquation(GraphScene, MovingCameraScene):
    CONFIG = {
        "x_min": -10, "x_max": 10, "x_axis_width": 20,
        "y_min": -5, "y_max": 5, "y_axis_height": 10,
        "axes_color": BLACK,
        "graph_origin": ORIGIN,
    }
    def setup(self):
        GraphScene.setup(self)
        MovingCameraScene.setup(self)
    
    def construct(self):
        self.setup_axes()
        self.play(FadeOut(self.axes))
        self.wait()
        nplane = VGroup(DoubleArrow(8 * LEFT, 15 * RIGHT), DoubleArrow(5.5 * DOWN, 8 * UP)).set_color(GREY)
        planegrid = VGroup()
        for i in range(-10, 15 + 1):
            planegrid.add(Line(8 * UP + i * RIGHT, 5 * DOWN + i * RIGHT))
        for i in range(-5, 8 + 1):
            planegrid.add(Line(i * UP + 15 * RIGHT, i * UP + 10 * LEFT))
        planegrid.set_color(GREY)
        planegrid.fade(0.875)
        spara = ParametricFunction(lambda t: np.array([t, 2 * t * np.sqrt(t) / 3, 0]), t_min = 0, t_max = 5, color = PURPLE)
        arealabel = TexMobject("A(x)", "=", "\\frac{2}{3}x^{3/2}")
        arealabel.move_to(3 * LEFT + UP + arealabel.get_center() - arealabel[1].get_center())
        arealabel.shift(UP)
        hparab = self.get_graph(
            lambda x: x ** 2,
            x_min = 0,
            x_max = 3,
            color = BLUE,
        ).flip(UR, about_point = ORIGIN)
        hparab.flip(RIGHT, about_point = ORIGIN)
        self.camera_frame.shift(UP + RIGHT)
        scub = self.get_graph(
            lambda x: 2 * np.sqrt(x ** 3) / 3,
            x_min = 0,
            x_max = 10,
            color = PURPLE,
        )
        auxcur = self.get_graph(
            lambda x: np.sqrt(1 + x),
            x_min = 0,
            x_max = 15,
            color = YELLOW,
        )
        self.add(nplane, planegrid, hparab, scub, auxcur)
        t = 3
        c = Circle(arc_center = t * RIGHT)
        pdot = Dot().move_to(t * RIGHT + 2 * np.sqrt(t ** 3) * UP / 3)
        plab = TexMobject("P").scale(0.75).add_background_rectangle().next_to(pdot, RIGHT, buff = 0.05)
        odot = Dot().move_to(t * RIGHT)
        olab = TexMobject("O").scale(0.75).add_background_rectangle().next_to(odot, DR, buff = 0.0)
        qdot = Dot().move_to(t * RIGHT + RIGHT)
        qlab = TexMobject("Q").scale(0.75).add_background_rectangle().next_to(qdot, DR, buff = 0.0)
        pddot = Dot().move_to(t * RIGHT + np.sqrt(t) * DOWN)
        pdlab = TexMobject("P'").scale(0.75).add_background_rectangle().next_to(pddot, UR, buff = 0.0)
        vline = Line(10 * UP, 10 * DOWN, color = GREEN).shift(t * RIGHT)
        tdot = Dot().move_to(t * RIGHT + np.sqrt(1 + t) * UP)
        tlab = TexMobject("T").scale(0.75).add_background_rectangle().next_to(tdot, DR, buff = 0.0)
        radline = Line(ORIGIN, UP, color = ORANGE).rotate_about_origin(np.arctan(np.sqrt(t))).shift(t * RIGHT)
        def lineextend(s0, e0, color = RED):
            s, e = s0, e0
            temp = Line(s0, e0)
            leng = temp.get_length()
            return Line(s0 + 100 * (e0 - s0) / leng, s0 - 100 * (e0 - s0) / leng, color = color)
        qpdline = lineextend(qdot.get_center(), pddot.get_center()).set_color(ORANGE)
        tddot = Dot().move_to(radline.get_end())
        tdlab = TexMobject("T'").scale(0.75).add_background_rectangle().next_to(tddot, UL, buff = 0.0)
        ctang = lineextend(tdot.get_center(), tddot.get_center()).set_color(ORANGE)
        rangle = Angle(tdot.get_center(), tddot.get_center(), odot.get_center(), radius = 0.25, color = YELLOW)
        self.play(*[Write(obj) for obj in [c, vline, tdot, tlab, pdot, plab, odot, olab, qdot, qlab, pddot, pdlab]])
        self.wait()
        self.play(*[Write(obj) for obj in [qpdline, ctang, radline, rangle, tddot, tdlab]])
        self.wait()
        eqn = TexMobject("{TO", "\\over", "OT'}", "=", "{P'Q", "\\over", "QO}").move_to(3 * UL)
        self.play(Write(eqn))
        self.wait()
        fone = TexMobject("1}").move_to(eqn[2].get_center())
        sone = TexMobject("1}").move_to(eqn[-1].get_center())
        self.play(
            Transform(eqn[2], fone),
            Transform(eqn[-1], sone),
        )
        self.wait()
        self.play(
            *[FadeOut(obj) for obj in [radline, rangle, tddot, tdlab, ctang, eqn[1], eqn[2], eqn[-1], eqn[-2]]],
            eqn[0].next_to, eqn[3], {"direction": LEFT},
            eqn[-3].next_to, eqn[3], {"direction": RIGHT},
        )
        self.play(
            eqn[0].shift, LEFT,
            eqn[3].shift, LEFT,
            eqn[-3].shift, LEFT,
        )
        self.wait()
        pyth = TexMobject("=", "\\sqrt{OQ^2+OP'^2}")
        pyth.move_to(eqn[3].get_center() + pyth.get_center() - pyth[0].get_center())
        pyth.shift(1 * DOWN)
        #pyth.next_to(eqn[3].get_center(), DOWN)
        self.play(Write(pyth))
        self.wait()
        parabeqn = TexMobject("y=", "-\\sqrt{x}").move_to(6 * RIGHT + 1.5 * DOWN)
        parabeqn[-1].set_color(BLUE)
        auxeq = TexMobject("=", "\\sqrt{1+(-\\sqrt{x})^2}")
        auxeq.move_to(pyth[0].get_center() + auxeq.get_center() - auxeq[0].get_center())
        auxeq.shift(1 * DOWN)
        self.play(
            Write(parabeqn),
            Write(auxeq)
        )
        self.wait()
        fineqn = TexMobject("=", "\\sqrt{1+x}")
        fineqn.move_to(auxeq[0].get_center() + fineqn.get_center() - fineqn[0].get_center())
        self.play(
            FadeOutAndShiftDown(auxeq),
            FadeIn(fineqn)
        )
        self.wait()
        self.play(
            *[FadeOut(obj) for obj in [eqn[3], eqn[-3], pyth, qdot, qlab, pddot, pdlab, qpdline, c]],
            fineqn.next_to, eqn[0], {"direction": RIGHT}
        )
        self.wait()
        pcpy = hparab.copy()
        self.play(pcpy.flip, RIGHT, {"about_point": ORIGIN})
        for _ in range(2):
            self.play(pcpy.shift, LEFT)
            self.play(pcpy.shift, RIGHT)
        self.play(
            FadeOut(pcpy),
            FadeOut(eqn[0]),
            FadeOut(fineqn),
        )
        pcpy.shift(LEFT)
        self.wait()
        def areacreator(minx = 0, maxx = 2, colorpref = YELLOW):
            tarea = VMobject(fill_color = colorpref, fill_opacity = 0.5, stroke_width = 0)
            areapts = [np.array([t, np.sqrt(1 + t), 0]) for t in np.arange(minx, maxx + 0.05, 0.05)]
            areapts = [np.array([minx, 0, 0])] + areapts + [np.array([maxx, 0, 0])]
            tarea.set_points_as_corners(areapts)
            return tarea
        #farea = areacreator(minx = -1, maxx = t)
        oarea = areacreator(minx = -1, maxx = 0)
        rarea = areacreator(minx = 0, maxx = t)
        self.play(Write(rarea))
        self.wait()
        self.play(Write(oarea))
        self.wait()
        self.play(FadeOut(oarea))
        self.wait()
        arclen = TexMobject("\\text{length}", "=", "\\frac{2}{3}(1+x)^{3/2}", "-", "\\frac{2}{3}").scale(0.875).add_background_rectangle()
        arclen.move_to(UP + 3 * LEFT)
        self.play(Write(arclen))
        self.wait(5)

if __name__ == "__main__":
    module_name = os.path.basename(__file__)
    #command_A = "manim -p  --video_dir ~/Downloads/  "
    clear_cmd = "cls"
    command_A = "python -m manim "
    #command_B = module_name + " " + "AuxEquation" + " -pl -n 2"
    #command_B = module_name + " " + "ProblemIntroduction" + " -ps -n 16,17"
    #command_B = module_name + " " + "IntroductionScene" + " -p",
    command_B = module_name + " -p"
    #command_B = module_name + " -pl"
    #command_B = module_name + " -a"
    #command_B = module_name + " -al"
    os.system(clear_cmd)
    os.system(command_A + command_B)