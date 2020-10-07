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

class EndScreen(GraphScene, MovingCameraScene):
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
        test = ParametricFunction(lambda t: np.array([t - np.sin(t) - PI, 1 - np.cos(t), 0]), t_min = -10, t_max = 10, color = BLUE)
        self.play(Write(test))
        self.play(
            self.camera_frame.move_to, 1 * UP,
            self.camera_frame.scale, 0.625,
        )
        cir = Circle().shift(UP)
        self.play(Write(cir))
        def cycloid_tangents(x_min = 0.0, x_max = 2, start_color = BLUE, end_color = GREEN, dx = 0.1):
            tanlist = VGroup()
            x_range = np.arange(x_min + dx, x_max + dx, dx)
            colors = color_gradient([start_color, end_color], len(x_range))
            cyc = lambda t: np.array([t - np.sin(t), 1 - np.cos(t), 0])
            for x, color in zip(x_range, colors):
                t = np.arccos(1 - x)
                if (2 - (1 - np.cos(t))) == 0:
                    tl = (PI - temp.get_start()[0]) * LEFT
                else:
                    tl = dx * (t - (t - np.sin(t))) / (2 - (1 - np.cos(t))) * LEFT + dx * DOWN
                temp = Line(cyc(t), cyc(t) + tl)
                temp.set_color(color)
                tanlist.add(temp)
                #tanlist.add(Dot(radius = 0.05, point = temp.get_start()).set_color(color))
            return tanlist
        def cycloid_tangent_list(max_dx = 1.0, n_iterations = 8):
            return [cycloid_tangents(dx = float(max_dx) / 2 ** n).shift(PI * LEFT) for n in range(1, n_iterations)]
        tanlist = cycloid_tangent_list()
        rtanlist = VGroup()
        for obj in tanlist:
            rtanlist.add(obj.copy().flip(about_point = ORIGIN))
        newtans = tanlist[0]
        rnewtans = rtanlist[0]
        self.play(
            ShowCreation(newtans),
            ShowCreation(rnewtans),
            #Transform(tlist, tlist2, run_time = 2, lag_ratio = 0.25),
            test.fade, 0.875,
            #self.camera_frame.scale, 0.05,
            #self.camera_frame.move_to, test(PI / 2),
        )
        self.wait()
        for i in range(1, len(tanlist)):
            self.play(
                Transform(newtans, tanlist[i], lag_ratio = 0.25),
                Transform(rnewtans, rtanlist[i], lag_ratio = 0.25)
            )
            self.wait()
        self.play(FadeOut(newtans), FadeIn(tanlist[2]), FadeOut(rnewtans), FadeIn(rtanlist[2]))
        self.wait(5)

class CycloidTangent(GraphScene, MovingCameraScene):
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
        planegrid = VGroup()
        for i in range(-10, 10 + 1):
            planegrid.add(Line(5 * UP + i * RIGHT, 5 * DOWN + i * RIGHT))
        for i in range(-5, 5 + 1):
            planegrid.add(Line(i * UP + 10 * RIGHT, i * UP + 10 * LEFT))
        planegrid.set_color(GREY)
        planegrid.fade(0.875)
        self.play(Write(nplane))
        self.play(ShowCreation(planegrid), run_time = 3)
        self.wait()
        cyc = lambda t: np.array([t - np.sin(t), 1 - np.cos(t), 0])
        cycloid = ParametricFunction(cyc, t_min = -10, t_max = 10, color = BLUE)
        cycloid.shift(PI * LEFT)
        self.play(Write(cycloid))
        self.wait()
        self.play(
            self.camera_frame.scale, 0.625,
            self.camera_frame.shift, UP,
        )
        gengrp = VGroup(
            Circle(),
            Dot(radius = 0.625 * DEFAULT_DOT_RADIUS, point = DOWN),
            Line(ORIGIN, DOWN)
        )
        gengrp.set_color(PURPLE)
        gengrp.shift(PI * LEFT + UP)
        self.play(Write(gengrp))
        pos = ValueTracker(-PI)
        def updfunc(obj):
            c, d, l = obj
            p = pos.get_value()
            tgp = VGroup(Circle(), Dot(radius = 0.625 * DEFAULT_DOT_RADIUS, point = DOWN), Line(ORIGIN, DOWN))
            tgp.set_color(PURPLE)
            tgp.rotate_about_origin(-p + PI)
            tgp.shift(UP + p * RIGHT)
            obj.become(tgp)
        tandot = Dot(point = cyc(PI / 2 + PI / 6) + PI * LEFT, radius = 0.625 * DEFAULT_DOT_RADIUS, stroke_width = 0.5, color = PURPLE)
        self.play(Write(tandot))
        self.wait()
        self.play(
            pos.set_value, -PI / 3,
            UpdateFromFunc(gengrp, updfunc), rate_func = smooth, run_time = 5)
        self.wait()
        icenter = Dot(color = YELLOW, radius = 0.625 * DEFAULT_DOT_RADIUS, point = gengrp[0].get_center() + DOWN)
        self.play(Write(icenter), FadeOut(gengrp[-1]))
        self.wait()
        theta = ValueTracker(225 * DEGREES)
        def arrowgrpcreator(t):
            tline = Arrow(DOWN, np.cos(t) * RIGHT + np.sin(t) * UP)
            rotgrp = VGroup(
                Dot(point = np.cos(t) * RIGHT + np.sin(t) * UP, radius = 0.625 * DEFAULT_DOT_RADIUS),
            )
            temp = Line(DOWN, np.cos(t) * RIGHT + np.sin(t) * UP)
            temp.rotate(PI / 2, about_point = temp.get_end())
            rotgrp.add(Angle(temp.get_start(), np.cos(t) * RIGHT + np.sin(t) * UP, DOWN, radius = 0.15, color = YELLOW))
            rotgrp.add(Arrow(temp.get_end(), temp.get_end() + (temp.get_start() - temp.get_end()) / temp.get_length()))
            rotgrp.add(tline)
            rotgrp.shift(icenter.get_center() + UP)
            return rotgrp
        def arrowgrpupdater(obj):
            thet = theta.get_value()
            obj.become(arrowgrpcreator(thet))
        rgrp = arrowgrpcreator(theta.get_value())
        self.play(Write(rgrp))
        self.play(theta.increment_value, -60 * DEGREES, UpdateFromFunc(rgrp, arrowgrpupdater), run_time = 2)
        self.play(theta.increment_value, -60 * DEGREES, UpdateFromFunc(rgrp, arrowgrpupdater), run_time = 2)
        self.play(theta.increment_value, -90 * DEGREES, UpdateFromFunc(rgrp, arrowgrpupdater), run_time = 2)
        self.wait()
        self.play(theta.set_value, 150 * DEGREES, UpdateFromFunc(rgrp, arrowgrpupdater), run_time = 2)
        self.wait()
        #self.play(FadeOut(rgrp))
        connline, tanline = Line(rgrp[-1].get_start(), rgrp[-1].get_end()), Line(rgrp[-2].get_start() + 10 * (rgrp[-2].get_end() - rgrp[-2].get_start()), rgrp[-2].get_start() - 10 * (rgrp[-2].get_end() - rgrp[-2].get_start()))
        self.play(ShowCreation(connline))
        self.wait()
        self.play(ShowCreation(tanline))
        self.wait()
        self.play(FadeOut(rgrp))
        self.wait()
        dialine = DashedLine(icenter.get_center(), icenter.get_center() + 2 * UP, color = GREY)
        tanseg = Line(tandot.get_center(), icenter.get_center() + 2 * UP)
        diadot = Dot(radius = 0.625 * DEFAULT_DOT_RADIUS, point = dialine.get_end())
        self.play(
            connline.fade, 0.5,
            tanline.fade, 0.5,
            ShowCreation(dialine),
            Write(diadot)
        )
        self.wait()
        self.play(Write(tanseg))
        self.wait(5)

class WrenFirstPart(GraphScene, MovingCameraScene):
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
        planegrid = VGroup()
        for i in range(-10, 10 + 1):
            planegrid.add(Line(5 * UP + i * RIGHT, 5 * DOWN + i * RIGHT))
        for i in range(-5, 5 + 1):
            planegrid.add(Line(i * UP + 10 * RIGHT, i * UP + 10 * LEFT))
        planegrid.set_color(GREY)
        planegrid.fade(0.875)
        cyc = lambda t: np.array([t - np.sin(t), 1 - np.cos(t), 0])
        cycloid = ParametricFunction(cyc, t_min = -10, t_max = 10, color = BLUE)
        cycloid.shift(PI * LEFT)
        self.wait()
        self.camera_frame.scale(0.625)
        self.camera_frame.shift(UP)
        self.add(
            nplane,
            planegrid,
            cycloid,
        )
        self.wait()
        def cycloid_tangents(x_min = 0.0, x_max = 2, start_color = BLUE, end_color = GREEN, dx = 0.1):
            tanlist = VGroup()
            x_range = np.arange(x_min + dx, x_max + dx, dx)
            colors = color_gradient([start_color, end_color], len(x_range))
            cyc = lambda t: np.array([t - np.sin(t), 1 - np.cos(t), 0])
            for x, color in zip(x_range, colors):
                t = np.arccos(1 - x)
                if (2 - (1 - np.cos(t))) == 0:
                    tl = (PI - temp.get_start()[0]) * LEFT
                else:
                    tl = dx * (t - (t - np.sin(t))) / (2 - (1 - np.cos(t))) * LEFT + dx * DOWN
                temp = Line(cyc(t), cyc(t) + tl)
                temp.set_color(color)
                tanlist.add(temp)
                #tanlist.add(Dot(radius = 0.05, point = temp.get_start()).set_color(color))
            return tanlist
        def cycloid_tangent_list(max_dx = 1.0, n_iterations = 7):
            return [cycloid_tangents(dx = float(max_dx) / 2 ** n).shift(PI * LEFT) for n in range(1, n_iterations)]
        tanlist = cycloid_tangents(dx = 0.25).shift(PI * LEFT)
        n = 3
        tseg = tanlist[n].copy()
        sdot = Dot(tanlist[n].get_start(), radius = 0.625 * DEFAULT_DOT_RADIUS)
        edot = Dot(tanlist[n - 1].get_start(), radius = 0.625 * DEFAULT_DOT_RADIUS)
        st, et = sdot.get_center()[1], edot.get_center()[1]
        st, et = np.arccos(1 - st), np.arccos(1 - et)
        cycseg = ParametricFunction(cyc, t_min = st, t_max = et, color = BLUE)
        cycseg.shift(PI * LEFT)
        self.play(Write(sdot), Write(edot))
        self.wait()
        self.play(
            ShowCreation(tanlist),
            cycloid.fade, 0.875
        )
        self.wait()
        self.play(
            FadeOut(tanlist),
            FadeIn(tseg),
            FadeIn(cycseg),
            FadeOut(sdot), FadeOut(edot),
            #self.camera_frame.move_to, edot.get_center(),
            self.camera_frame.scale, 0.75
        )
        self.wait()
        templen = tanlist[n].get_length()
        templen = templen *  (2 - tanlist[n].get_start()[1]) / (tanlist[n].get_start()[1] - tanlist[n].get_end()[1])
        stan = Line(tanlist[n].get_start(), tanlist[n].get_start() - templen * (tanlist[n].get_end() - tanlist[n].get_start()) / tanlist[n].get_length())
        templen = tanlist[n - 1].get_length()
        templen = templen *  (2 - tanlist[n - 1].get_start()[1]) / (tanlist[n - 1].get_start()[1] - tanlist[n - 1].get_end()[1])
        etan = Line(tanlist[n - 1].get_start(), tanlist[n - 1].get_start() - templen * (tanlist[n - 1].get_end() - tanlist[n - 1].get_start()) / tanlist[n - 1].get_length())
        stan.set_color(YELLOW)
        etan.set_color(YELLOW)
        gcir = Circle().shift(UP)
        ecir, scir = Circle(arc_center = etan.get_end() + DOWN), Circle(arc_center = stan.get_end() + DOWN)
        egp, sgp = VGroup(ecir, DashedLine(etan.get_end(), etan.get_end() + 2 * DOWN).set_color(GREY)), VGroup(scir, DashedLine(stan.get_end(), stan.get_end() + 2 * DOWN).set_color(GREY))
        egp.fade(0.5)
        sgp.fade(0.5)
        self.play(ShowCreation(etan), ShowCreation(egp))
        self.wait()
        self.play(
            etan.shift, etan.get_end()[0] * LEFT,
            FadeOutAndShift(egp, etan.get_end()[0] * LEFT))
        self.wait()
        self.play(Write(stan), ShowCreation(gcir), ShowCreation(sgp))
        self.wait()
        self.play(
            stan.shift, stan.get_end()[0] * LEFT,
            tseg.shift, stan.get_end()[0] * LEFT,
            FadeOutAndShift(sgp, stan.get_end()[0] * LEFT))
        self.wait()
        #ecir, scirc = Circle(arc_center = 2 * UP, radius = etan.get_length, color = PURPLE), Circle(arc_center = 2 * UP, radius = stan.get_length, color = PURPLE)
        ecir, scir = Circle(arc_center = 2 * UP, radius = etan.get_length(), color = PURPLE), Circle(arc_center = 2 * UP, radius = stan.get_length(), color = PURPLE)
        self.play(Write(ecir), Write(scir))
        self.wait()
        eseg = Line(etan.get_end() + stan.get_length() * (etan.get_start() - etan.get_end()) / etan.get_length(), etan.get_start()).set_color(YELLOW)
        self.add(eseg)
        self.play(
            etan.fade, 0.75,
            stan.fade, 0.75,
        )
        purptan = Line(ORIGIN, 5 * etan.get_start()).set_color(PURPLE)
        redtan = Line(UP, etan.get_start())
        redtan.rotate(PI / 2, about_point = etan.get_start())
        redtan = Line(redtan.get_start() + 5 * (redtan.get_end() - redtan.get_start()), redtan.get_start() - 5 * (redtan.get_end() - redtan.get_start())).set_color(RED)
        self.play(
            Write(purptan),
            ecir.fade, 0.5)
        self.wait()
        self.play(
            Write(redtan),
            self.camera_frame.move_to, etan.get_start(),
            self.camera_frame.scale, 0.75)
        self.wait(5)

class AlmostCycloid(GraphScene, MovingCameraScene):
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
        planegrid = VGroup()
        for i in range(-10, 10 + 1):
            planegrid.add(Line(5 * UP + i * RIGHT, 5 * DOWN + i * RIGHT))
        for i in range(-5, 5 + 1):
            planegrid.add(Line(i * UP + 10 * RIGHT, i * UP + 10 * LEFT))
        planegrid.set_color(GREY)
        planegrid.fade(0.875)
        #self.play(Write(nplane))
        #self.play(ShowCreation(planegrid), run_time = 3)
        #self.wait()
        def cyc_generate(rtime = 10, colorpref = GREEN):
            t = ValueTracker(0)
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
            self.play(t.set_value, 2 * PI, run_time = rtime, rate_func = linear)
            cg.clear_updaters()
            pathgp.clear_updaters()
            self.play(FadeOut(cg), FadeOut(pathd))
            self.remove(path)
            return path
        def almost_cyc_generate(n, rtime = 10, colorpref = GREEN):
            t = ValueTracker(0)
            angspan = 180 * DEGREES / n
            c = RegularPolygon(n, start_angle = -90 * DEGREES - angspan)
            l = Line(ORIGIN, c.get_start())
            cg = VGroup(c, l)
            sidelen = Line(c.points[0], c.points[3]).get_length()
            shiftdist = -c.points[0]
            cg.shift(shiftdist)
            path = VMobject(color = colorpref)
            pathd = Dot(cg[-1].get_end(), radius = DEFAULT_SMALL_DOT_RADIUS)
            path.set_points_as_corners([pathd.get_center(), pathd.get_center() + 0.001 * UP])
            pathgp = VGroup(path, pathd)
            def movcircle(obj):
                tc, tl = obj
                tempc = RegularPolygon(n, start_angle = -90 * DEGREES - angspan)
                templ = Line(ORIGIN, tempc.get_start())
                tempg = VGroup(tempc, templ)
                j = t.get_value()
                div = 360 * DEGREES / n
                quo, rem = j // div, j % div
                tempg.rotate_about_origin(-quo * 360 * DEGREES / n)
                tempg.shift(shiftdist)
                tempg.shift(quo * sidelen * RIGHT)
                tempg.rotate(-rem, about_point = (quo + 1) * sidelen * RIGHT)
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
            self.play(t.set_value, 2 * PI * (1 - 1 / n), run_time = rtime, rate_func = linear)
            cg.clear_updaters()
            pathgp.clear_updaters()
            self.play(FadeOut(cg), FadeOut(pathd))
            self.remove(path)
            return path
        self.play(
            self.camera_frame.scale, 0.75,
            self.camera_frame.shift, PI * RIGHT + UP,
        )
        baseline = Arrow(10 * LEFT, 10 * RIGHT).set_color(GREY)
        self.play(Write(baseline))
        self.wait()
        titles = VGroup(
            TexMobject("\\underline{\\text{Cycloid of a rolling Circle}}"),
            TexMobject("\\underline{\\text{(almost) Cycloid of a rolling Square}}"),
            TexMobject("\\underline{\\text{(almost) Cycloid of a rolling Pentagon}}"),
            TexMobject("\\underline{\\text{(almost) Cycloid of a rolling 10-gon}}"),
            TexMobject("\\underline{\\text{(almost) Cycloid of a rolling 15-gon}}"),
            TexMobject("\\underline{\\text{(almost) Cycloid of a rolling 20-gon}}"),
            TexMobject("\\underline{\\text{(almost) Cycloid of a rolling Hexagon}}"),
        )
        for obj in titles:
            obj.scale(0.875)
        titles.move_to(PI * RIGHT + 3.5 * UP)
        txt = titles[0].copy()
        self.play(Write(txt))
        self.wait()
        curves = VGroup()
        cyccurve = cyc_generate()
        curves.add(cyccurve)
        self.add(cyccurve)
        self.play(FadeOut(cyccurve))
        sidelist = [10000, 4, 5, 10, 15, 20, 6] 
        for i in [1, 2, 3, 4, 5, 6]:
            self.play(Transform(txt, titles[i].copy()))
            self.wait()
            temp = almost_cyc_generate(sidelist[i], colorpref = GREEN)
            curves.add(temp)
            self.add(temp)
            self.play(FadeOut(temp))
            self.wait()
        self.wait(5)

class HexCycloid(GraphScene, MovingCameraScene):
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
        planegrid = VGroup()
        for i in range(-10, 10 + 1):
            planegrid.add(Line(5 * UP + i * RIGHT, 5 * DOWN + i * RIGHT))
        for i in range(-5, 5 + 1):
            planegrid.add(Line(i * UP + 10 * RIGHT, i * UP + 10 * LEFT))
        planegrid.set_color(GREY)
        planegrid.fade(0.875)
        def cyc_generate(rtime = 10, colorpref = GREEN):
            t = ValueTracker(0)
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
            self.play(t.set_value, 2 * PI, run_time = rtime, rate_func = linear)
            cg.clear_updaters()
            pathgp.clear_updaters()
            self.play(FadeOut(cg), FadeOut(pathd))
            self.remove(path)
            return path
        def almost_cyc_generate(n, rtime = 10, colorpref = GREEN):
            t = ValueTracker(0)
            angspan = 180 * DEGREES / n
            c = RegularPolygon(n, start_angle = -90 * DEGREES - angspan)
            l = Line(ORIGIN, c.get_start())
            cg = VGroup(c, l)
            sidelen = Line(c.points[0], c.points[3]).get_length()
            shiftdist = -c.points[0]
            cg.shift(shiftdist)
            path = VMobject(color = colorpref)
            pathd = Dot(cg[-1].get_end(), radius = DEFAULT_SMALL_DOT_RADIUS)
            path.set_points_as_corners([pathd.get_center(), pathd.get_center() + 0.001 * UP])
            pathgp = VGroup(path, pathd)
            def movcircle(obj):
                tc, tl = obj
                tempc = RegularPolygon(n, start_angle = -90 * DEGREES - angspan)
                templ = Line(ORIGIN, tempc.get_start())
                tempg = VGroup(tempc, templ)
                j = t.get_value()
                div = 360 * DEGREES / n
                quo, rem = j // div, j % div
                tempg.rotate_about_origin(-quo * 360 * DEGREES / n)
                tempg.shift(shiftdist)
                tempg.shift(quo * sidelen * RIGHT)
                tempg.rotate(-rem, about_point = (quo + 1) * sidelen * RIGHT)
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
            self.play(t.set_value, 2 * PI * (1 - 1 / n), run_time = rtime, rate_func = linear)
            cg.clear_updaters()
            pathgp.clear_updaters()
            self.play(FadeOut(cg), FadeOut(pathd))
            self.remove(path)
            return path
        self.play(
            self.camera_frame.scale, 0.75,
            self.camera_frame.shift, PI * RIGHT + UP,
        )
        baseline = Arrow(10 * LEFT, 10 * RIGHT).set_color(GREY)
        self.add(baseline)
        curv = almost_cyc_generate(6)
        ccurv = curv.copy()
        self.add(curv)
        poly = RegularPolygon(6, start_angle = -90 * DEGREES - 180 * DEGREES / 6)
        sidelen = Line(poly.points[0], poly.points[3]).get_length()
        poly.shift(-poly.points[0])
        self.wait()
        self.play(ShowCreation(poly))
        self.wait()
        self.play(
            self.camera_frame.shift, sidelen * RIGHT,
            poly.shift, 7 * sidelen * RIGHT,
        )
        t = RegularPolygon(6, start_angle = -90 * DEGREES - 180 * DEGREES / 6)
        t.shift(-t.points[0])
        lgp = VGroup()
        dgp = VGroup()
        for i in range(1, 6):
            tg = VGroup()
            temp = Line(t.points[0], t.points[3])
            temp.rotate((i - 1) * 360 * DEGREES / 6, about_point = t.get_center())
            temp = Line(t.points[0], temp.get_end())
            tempp = temp.copy()
            tempp.set_color(YELLOW)
            tempp.shift(7 * sidelen * RIGHT)
            tempp.rotate(PI, axis = UP, about_point = 7.5 * sidelen * RIGHT)
            dgp.add(tempp)
            for j in range(1, i):
                temp.rotate(-2 * PI / 6, about_point = j * sidelen * RIGHT)
            rtemp = temp.copy()
            rtemp.rotate(-2 * PI / 6, about_point = temp.get_end())
            tg.add(temp)
            tg.add(rtemp)
            temp.set_color(YELLOW)
            rtemp.set_color(PURPLE)
            self.play(Write(tg), Write(tempp))
            self.wait()
            lgp.add(tg)
        anggp = VGroup()
        for obj in lgp:
            anggp.add(Angle(obj[0].get_start(), obj[0].get_end(), obj[1].get_start(), radius = 0.375))
        anggp.add(Angle(dgp[0].get_start() + RIGHT, dgp[0].get_start(), dgp[-1].get_end(), radius = 0.375))
        hexangle = TexMobject("2\\pi/6").scale(0.5)
        hexangle.next_to(anggp[-1], RIGHT, buff = 0.1)
        self.play(Write(anggp), Write(hexangle))
        self.wait(2)
        ccircle = Circle().shift(4 * RIGHT + 3 * DOWN).set_color(GREEN)
        polyg = RegularPolygon(8).shift(4 * RIGHT + 3 * DOWN)
        self.play(
            self.camera_frame.shift, 4 * DOWN,
            self.camera_frame.scale, 0.625,
            Write(ccircle),
            Write(polyg)
        )
        self.wait()
        line1 = Line(ccircle.get_center(), ccircle.get_start()).set_color(PURPLE)
        line2 = line1.copy().rotate(2 * PI / 8, about_point = line1.get_start())
        ang = Angle(line1.get_end(), line1.get_start(), line2.get_end(), radius = 0.25)
        sidegp = VGroup(line1, line2, ang)
        scpy = sidegp.copy()
        self.add(scpy)
        thet = TexMobject("\\theta").scale(0.375).next_to(ang, RIGHT, buff = 0.07)
        self.play(Write(sidegp), Write(thet))
        self.wait()
        for i in range(3):
            self.play(scpy.rotate, -360 * DEGREES / 8, {"about_point": ccircle.get_center()})
        self.play(
            Rotate(scpy, angle = -360 * DEGREES / 8, about_point = ccircle.get_center()),
            VFadeOut(scpy),
            FadeOut(thet),
            FadeOut(sidegp)
        )
        self.wait()
        sixteengon = VMobject(color = BLUE)
        temp = line1.copy()
        rang = 360 * DEGREES / 12
        pts = []
        temp.rotate(-2 * rang, about_point = ccircle.get_center())
        for i in range(11):
            pts += [temp.get_end()]
            temp.rotate(rang, about_point = ccircle.get_center())
        sixteengon.set_points_as_corners(pts)
        dline = DashedLine(pts[-1], pts[0], color = GREY)
        self.play(FadeOut(polyg), Write(sixteengon), Write(dline))
        self.wait()
        kangle = VGroup(
            Line(ccircle.get_center(), pts[1], color = PURPLE),
            Line(ccircle.get_center(), pts[-2], color = PURPLE),
            Line(pts[1], pts[-2], color = YELLOW),
            Angle(pts[1], ccircle.get_center(), pts[-2], radius = 0.1875)
        )
        ksidelabel = TexMobject("2\\sin(k\\pi/n)").scale(0.375).next_to(kangle[2], DOWN, buff = 0.05)
        ksidelabel.add_background_rectangle()
        kanglelabel = TexMobject("2k\\pi/n").scale(0.375).next_to(ccircle.get_center(), UP, buff = 0.05)
        self.play(ShowCreation(kangle))
        self.wait()
        self.play(Write(kanglelabel))
        self.wait()
        self.play(Write(ksidelabel))
        self.wait()
        self.play(
            *[FadeOut(obj) for obj in [kangle, kanglelabel, ksidelabel, ccircle, dline, sixteengon]],
            self.camera_frame.shift, 3 * UP,
            self.camera_frame.scale, 1 / 0.625,
        )
        arlen = VGroup(
            TextMobject("Arc Length,"),
            TexMobject("L_6"),
            TexMobject("="),
            TexMobject("\\displaystyle \\sum_{k=0}^6 2\\sin\\left(\\frac{k\\pi}{6}\\right)\\cdot \\frac{2\\pi}{6}")
        )
        for obj in arlen:
            obj.scale(0.875)
            obj.add_background_rectangle()
        arlen.arrange(direction = RIGHT, buff = 0.25)
        arlen.shift(4 * RIGHT + 1.5 * DOWN)
        self.play(Write(arlen))
        self.wait()
        lneqn = TexMobject("\\displaystyle \\sum_{k=0}^n 2\\sin\\left(\\frac{k\\pi}{n}\\right)\\cdot \\frac{2\\pi}{n}").scale(0.875).move_to(arlen[-1].get_center())
        self.play(
            Transform(arlen[1], TexMobject("L_n").scale(0.875).move_to(arlen[1].get_center())),
            Transform(hexangle, TexMobject("2\\pi/n").scale(0.5).next_to(anggp[-1], RIGHT, buff = 0.1)),
            FadeOutAndShiftDown(arlen[-1]),
            FadeInFrom(lneqn, UP)
        )
        self.wait()
        integeqns = VGroup(
            TexMobject("L_\\infty"),
            TexMobject("=\\displaystyle \\lim_{n \\to \\infty}4\\pi \\sum_{k=0}^n \\sin\\left(\\frac{k\\pi}{n}\\right)\\cdot\\frac{1}{n}"),
            TexMobject("=\\displaystyle 4\\pi\\int\\limits_{0}^1 \\sin(\\pi x)\\,dx")
        )
        for obj in integeqns:
            obj.scale(0.875)
        integeqns[1].next_to(integeqns[0], RIGHT, buff = 0.25)
        integeqns[2].next_to(integeqns[0], RIGHT, buff = 0.25)
        integeqns.shift(1.5 * RIGHT + 4 * DOWN)
        self.play(
            self.camera_frame.shift, 3 * DOWN,
            Write(integeqns[:2])
        )
        self.wait()
        self.play(
            FadeOutAndShiftDown(integeqns[1]),
            FadeInFrom(integeqns[2], UP)
        )
        self.wait()
        cr = Cross(integeqns[2])
        self.play(ShowCreation(cr))
        self.wait()
        self.play(*[FadeOut(obj) for obj in [cr, integeqns[0], integeqns[2]]])
        self.wait()
        sumform = TexMobject("=\\displaystyle \\frac{4\\pi}{n}\\sum_{k=0}^n \\sin\\left(\\frac{k\\pi}{n}\\right)").scale(0.875)
        sumform.next_to(lneqn, DOWN, buff = 0.25)
        sumform.align_to(arlen[2], direction = LEFT)
        self.play(Write(sumform))
        self.wait()
        lagiden = VGroup(
            TexMobject("\\underline{\\text{Lagrange's Identity}}").scale(0.875).set_color(YELLOW),
            #TexMobject("\\displaystyle \\sum_{k=0}^n \\sin(k\\theta)=\\frac{\\sin\\left(\\frac{n\\theta}{2}\\right)\\sin\\left(\\frac{(n+1)\\theta}{2}\\right)}{\\sin\\left(\\frac{\\theta}{2}\\right)}\\cdot").scale(0.875)
            #TexMobject("\\displaystyle \\sum_{k=0}^n \\sin(k\\theta)=\\frac{\\sin(n\\theta / 2)\\sin(n\\theta / 2 + \\theta/2)}{\\sin(\\theta/2)}").scale(0.875)
            TexMobject("\\displaystyle \\sum_{k=0}^n \\sin(k\\theta)=\\frac{\\sin(\\frac{n\\theta}{2})\\sin(\\frac{n\\theta}{2}+\\frac{\\theta}{2})}{\\sin(\\frac{\\theta}{2})}").scale(0.875)
        )
        lagiden.arrange(direction = DOWN, buff = 0.25)
        #lagiden[1].align_to(lagiden[0], direction = LEFT)
        lagiden.shift(4 * RIGHT + 4 * DOWN)
        self.play(
            FadeOutAndShift(arlen[2], UP),
            FadeOutAndShift(lneqn, UP),
            sumform.next_to, arlen[1], {"direction": RIGHT, "buff": 0.25})
        self.wait()
        self.play(Write(lagiden))
        self.wait()
        finaleqn = VGroup(
            TexMobject("L_\\infty"),
            TexMobject("=\\displaystyle \\lim_{n \\to \\infty}\\frac{4\\pi}{n}\\sum_{k=0}^n\\sin\\left(\\frac{k\\pi}{n}\\right)"),
            TexMobject("=\\displaystyle \\lim_{n \\to \\infty}\\frac{4\\pi}{n}\\frac{\\sin(\\frac{\\pi}{2})\\sin(\\frac{\\pi}{2} + \\frac{\\pi}{2n})}{\\sin(\\frac{\\pi}{2n})}"),
            TexMobject("=\\displaystyle \\lim_{n \\to \\infty}\\frac{4\\pi}{n}\\frac{1}{\\sin(\\frac{\\pi}{2n})}"),
            TexMobject("=\\displaystyle \\lim_{n \\to \\infty}8\\frac{\\frac{\\pi}{2n}}{\\sin(\\frac{\\pi}{2n})}"),
            TexMobject("=\\displaystyle 8"),
        )
        for text in finaleqn:
            text.scale(0.875)
        for text in finaleqn[1:]:
            text.next_to(finaleqn[0], RIGHT, buff = 0.25)
        finaleqn.shift(1 * RIGHT + 6.5 * DOWN)
        self.play(
            Write(finaleqn[:2]),
            self.camera_frame.shift, 2.5 * DOWN,
        )
        self.wait()
        for i in range(2, 1 + 4):
            self.play(
                FadeOutAndShiftDown(finaleqn[i - 1]),
                FadeInFrom(finaleqn[i], UP),
            )
            self.wait()
        ultxt = TexMobject("\\displaystyle \\lim_{x \\to 0}\\frac{\\sin x}{x}=1").scale(0.875).set_color(YELLOW)
        ulbub = ThoughtBubble()
        ulbub.rotate(PI, axis = UP)
        #ulbub.pin_to(finaleqn[-1])
        #ulbub.pin_to(Dot().to_corner(DR))
        ulbub.next_to(finaleqn[-2], direction = RIGHT)
        ulbub.add_content(ultxt)
        ulbub.resize_to_content()
        ulgrp = VGroup(ulbub, ultxt)
        self.play(Write(ulgrp))
        self.wait()
        self.play(
            FadeOutAndShiftDown(finaleqn[-2]),
            FadeInFrom(finaleqn[-1], UP)
        )
        self.wait()
        self.play(
            *[FadeOut(obj) for obj in [ulgrp, finaleqn[0], finaleqn[-1], lagiden]],
            Transform(arlen[1], TexMobject("L_\\infty").scale(0.875).move_to(arlen[1].get_center())),
            Transform(sumform, TexMobject("=\\displaystyle \\lim_{n \\to \\infty}\\frac{4\\pi}{n}\\sum_{k=0}^n \\sin\\left(\\frac{k\\pi}{n}\\right)").scale(0.875).next_to(arlen[1], buff = 0.25)),
            self.camera_frame.shift, 5.5 * UP)
        self.wait()
        cyceqn = ParametricFunction(lambda t: np.array([t - np.sin(t), 1 - np.cos(t), 0]), t_min = 0, t_max = 2 * PI, color = GREEN)
        fcircle = Circle(color = BLUE).shift(poly.get_center()[0] * RIGHT + UP)
        tgp = VGroup(arlen[0], arlen[1], sumform)
        finalresult = TextMobject("Arc Length = 8 $\\times$ radius of the generating circle").scale(0.875).move_to(4 * RIGHT + 3.5 * UP)
        finalresult[0][0:9].set_color(GREEN)
        finalresult[0][-16:].set_color(BLUE)
        self.play(
            FadeOut(lgp),
            FadeOut(dgp),
            FadeOut(anggp),
            FadeOut(hexangle),
            Transform(curv, cyceqn),
            Transform(poly, fcircle),
            FadeOutAndShiftDown(tgp),
            Write(finalresult),
            self.camera_frame.shift, 1.5 * UP
        )
        self.wait()
        self.play(
            poly.move_to, PI * RIGHT + UP,
            finalresult.shift, 0.75 * LEFT,
            self.camera_frame.shift, 0.75 * LEFT
        )
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
        self.wait()
        excurve = ParametricFunction(lambda t: np.array([np.exp(0.2 * t) * np.cos(t), np.exp(0.2 * t) * np.sin(t), 0]), t_min = -25, t_max = 10, color = BLUE)
        self.play(ShowCreation(excurve), run_time = 10, rate_func = linear)
        self.wait()
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
            self.play(t.set_value, 3 * PI, run_time = rtime, rate_func = linear)
            cg.clear_updaters()
            pathgp.clear_updaters()
            self.play(FadeOut(cg), FadeOut(pathd))
            self.remove(path)
            return path
        cycurve = cyc_generate()
        self.add(cycurve)
        self.wait()
        self.play(FadeOut(excurve))
        self.wait()
        mathnames = VGroup(
            TextMobject("Blaise Pascal"),
            TextMobject("Gilles de Roberval"),
            TextMobject("Evangelista Torricelli"),
            TextMobject("Christiaan Huygens")
        )
        for text in mathnames:
            #text.scale(0.9375)
            text.add_background_rectangle()
        #mathnames.arrange_in_grid(2, 2)
        mathnames[0].move_to(4 * LEFT + DOWN)
        mathnames[1].move_to(4 * RIGHT + DOWN)
        mathnames[2].move_to(4 * LEFT + 3 * DOWN)
        mathnames[3].move_to(4 * RIGHT + 3 * DOWN)
        for text in mathnames:
            self.play(Write(text))
        wren = ImageMobject("Wren").scale(2.5)
        wren.shift(UP)
        wrenname = VGroup(
            TextMobject("Christopher Wren"),
            TextMobject("(1632 - 1723)")
        )
        wrenname.arrange(DOWN)
        wrenname.next_to(wren, DOWN)
        wrenname.add_background_rectangle()
        self.play(
            FadeOut(mathnames),
            Write(wrenname),
            FadeIn(wren)
        )
        self.wait(5)

if __name__ == "__main__":
    module_name = os.path.basename(__file__)
    #command_A = "manim -p  --video_dir ~/Downloads/  "
    clear_cmd = "cls"
    command_A = "python -m manim "
    #command_B = module_name + " " + "Intro" + " -pl"
    #command_B = module_name + " " + "ProblemIntroduction" + " -ps -n 16,17"
    #command_B = module_name + " " + "IntroductionScene" + " -p",
    command_B = module_name + " -p"
    #command_B = module_name + " -a"
    #command_B = module_name + " -al"
    os.system(clear_cmd)
    os.system(command_A + command_B)