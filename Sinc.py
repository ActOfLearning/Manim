from manimlib.imports import *
import numpy as np
from scipy import signal

'''class ProblemIntroOutro(GraphScene, MovingCameraScene):
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
        def convolvedgrph(obj1, obj2):
            #make sure the x-coordinates of both the objects are the same
            ob1 = [x[1] for x in obj1]
            ob2 = [x[1] for x in obj2]
            l = len(ob2)
            res = signal.convolve(ob1, ob2) / sum(ob2)
            l = (l - 1) // 2
            res = res[l:]
            res = res[:-l]
            res = [[x[0], y, 0] for x, y in zip(obj2, res)]
            return res
        eps = 1 / 100
        res = [np.array([x, 1 if -0.5 <= x <= 0.5 else 0, 0]) for x in np.arange(-10, 10 + eps, eps)]
        res_g = VMobject().set_points_as_corners(res)
        c_list = [BLUE, GREEN, TEAL, GOLD, RED, PURPLE, ORANGE, PINK, MAROON]
        self.play(
            self.camera_frame.scale, 1 / 3,
            #self.camera_frame.shift, UP
        )
        for i in range(2, 2):
            temp = [np.array([x, 1 if -1 / 2 <= x <= 1 / 2 else 0, 0]) for x in np.arange(-10, 10 + eps, eps)]
            res = convolvedgrph(res, temp)
            res_g = VMobject().set_points_smoothly(res)
            res_g.set_color(c_list[i % 9])
            self.play(Write(res_g))
        def get_idx(arr, x):
            return np.searchsorted(arr, x)
        def get_from_x_coordinate(obj, x):
            temp = obj.points
            indx = get_idx(temp[:, 0], x)
            return temp[indx]
        movsq = VGroup(
            Line(LEFT / 2, LEFT / 2 + UP),
            Line(LEFT / 2 + UP, RIGHT / 2 + UP),
            Line(RIGHT / 2 + UP, RIGHT / 2),
            DashedLine(ORIGIN, UP),
        )
        exppts = [np.array([x, 0 if x < 0 else np.exp(-x), 0]) for x in np.arange(-10, 10 + eps, eps)]
        exp_g = VMobject().set_points_as_corners(exppts)
        convexppts = convolvedgrph(exppts, res)
        convexp_g = VMobject().set_points_as_corners(convexppts)
        def area_creator(obj1, x1, x2):
            resobj = VMobject(fill_color = YELLOW, fill_opacity = 0.5)
            tpts = obj1.points.copy()
            idx1, idx2 = get_idx(tpts[:, 0], x1), get_idx(tpts[:, 0], x2)
            tpts = tpts[idx1: idx2]
            tpts = np.insert(tpts, 0, [[x1, 0, 0]], axis = 0)
            tpts = np.append(tpts, [[x2, 0, 0], [x1, 0, 0]], axis = 0)
            resobj.set_points_as_corners(tpts)
            return resobj
        movsq.shift(2 * LEFT)
        exparea = area_creator(exp_g, movsq.get_critical_point(DOWN)[0] - 0.5, movsq.get_critical_point(DOWN)[0] + 0.5)
        exppath = VMobject(color = RED)
        expdot = Dot(movsq[-1].get_start(), color = RED, radius = 1 / 50)
        exppath.set_points_as_corners([expdot.get_center(), expdot.get_center() + 0.001 * UP])
        exppathgrp = VGroup(exppath, expdot, convexp_g.fade(1))
        exp_gp = VGroup(movsq, exparea, exp_g)
        self.play(Write(exp_gp), Write(exppathgrp))
        def conv_updater(obj):
            tm, ta, tc = obj
            curr = tm.get_critical_point(DOWN)[0]
            ta.become(area_creator(tc, curr - 0.5, curr + 0.5))
        def path_updater(obj):
            pat, pd, pc = obj
            curr = get_from_x_coordinate(pc, movsq[-1].get_start()[0])
            pd.move_to(curr)
            old_path = pat.copy()
            old_path.append_vectorized_mobject(Line(old_path.points[-1], curr))
            #old_path.make_smooth()
            pat.become(old_path)
        exp_gp.add_updater(conv_updater)
        exppathgrp.add_updater(path_updater)
        self.add(exp_gp, exppathgrp)
        self.play(movsq.shift, 4 * RIGHT, run_time = 10, rate_func = linear)
        exp_gp.clear_updaters()
        exppathgrp.clear_updaters()
        self.play(Write(convexp_g))
        self.wait(5)

class SawToothConvolution(GraphScene, MovingCameraScene):
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
        def convolvedgrph(obj1, obj2):
            #make sure the x-coordinates of both the objects are the same
            ob1 = [x[1] for x in obj1]
            ob2 = [x[1] for x in obj2]
            l = len(ob2)
            res = signal.convolve(ob1, ob2) / sum(ob2)
            l = (l - 1) // 2
            res = res[l:]
            res = res[:-l]
            res = [[x[0], y, 0] for x, y in zip(obj2, res)]
            return res
        eps = 1 / 50
        self.play(self.camera_frame.scale, 1 / 2)
        unitbox = [np.array([x, 1 if -0.5 <= x <= 0.5 else 0, 0]) for x in np.arange(-5, 5 + eps, eps)]
        unitbox_g = VMobject().set_points_as_corners(unitbox)
        #sawtooth = [np.array([x, x - np.floor(x), 0]) for x in np.arange(-10, 10 + eps, eps)]
        #sawtooth = [np.array([x, 0 if x == 0 else np.ceil(abs(x)) * ((abs(x) % 1) % (1 / np.ceil(abs(x)))), 0]) for x in np.arange(-10, 10 + eps, eps)]
        #sawtooth = [np.array([x, (np.sin(x) + np.sin(2 * x)) / 2, 0]) for x in np.arange(-7, 7 + eps, eps)]
        sawtooth = [np.array([x, np.exp(-abs(x)), 0]) for x in np.arange(-5, 5 + eps, eps)]
        sawtooth_g = VMobject().set_points_as_corners(sawtooth)
        def get_idx(arr, x):
            return np.searchsorted(arr, x)
        def get_from_x_coordinate(obj, x):
            temp = obj.points
            indx = get_idx(temp[:, 0], x)
            return temp[indx]
        movsq = VGroup(
            Line(LEFT / 2, LEFT / 2 + UP),
            Line(LEFT / 2 + UP, RIGHT / 2 + UP),
            Line(RIGHT / 2 + UP, RIGHT / 2),
            DashedLine(ORIGIN, UP),
        )
        convpts = convolvedgrph(sawtooth, unitbox)
        convpts_g = VMobject().set_points_as_corners(convpts)
        def area_creator(obj1, x1, x2):
            resobj = VMobject(fill_color = YELLOW, fill_opacity = 0.5)
            tpts = obj1.points.copy()
            idx1, idx2 = get_idx(tpts[:, 0], x1), get_idx(tpts[:, 0], x2)
            tpts = tpts[idx1: idx2]
            tpts = np.insert(tpts, 0, [[x1, 0, 0]], axis = 0)
            tpts = np.append(tpts, [[x2, 0, 0], [x1, 0, 0]], axis = 0)
            resobj.set_points_as_corners(tpts)
            return resobj
        movsq.shift(3.5 * LEFT)
        saw_area = area_creator(sawtooth_g, movsq.get_critical_point(DOWN)[0] - 0.5, movsq.get_critical_point(DOWN)[0] + 0.5)
        saw_path = VMobject(color = RED)
        saw_dot = Dot(movsq[-1].get_start(), color = RED, radius = 1 / 50)
        saw_path.set_points_as_corners([saw_dot.get_center(), saw_dot.get_center() + 0.001 * UP])
        saw_path_grp = VGroup(saw_path, saw_dot, convpts_g.fade(1))
        saw_grp = VGroup(movsq, saw_area, sawtooth_g)
        self.play(Write(saw_grp), Write(saw_path_grp))
        def conv_updater(obj):
            tm, ta, tc = obj
            curr = tm.get_critical_point(DOWN)[0]
            ta.become(area_creator(tc, curr - 0.5, curr + 0.5))
        def path_updater(obj):
            pat, pd, pc = obj
            curr = get_from_x_coordinate(pc, movsq[-1].get_start()[0])
            pd.move_to(curr)
            old_path = pat.copy()
            old_path.append_vectorized_mobject(Line(old_path.points[-1], curr))
            pat.become(old_path)
        saw_grp.add_updater(conv_updater)
        saw_path_grp.add_updater(path_updater)
        self.add(saw_grp, saw_path_grp)
        self.play(movsq.shift, 7 * RIGHT, run_time = 20, rate_func = linear)
        saw_grp.clear_updaters()
        saw_path_grp.clear_updaters()
        self.wait(5)
'''

class SincWaveFourier(GraphScene, MovingCameraScene):
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
        sinc = self.get_graph(lambda z: 1 if z == 0 else np.sin(z) / z)
        scaledsinc = self.get_graph(lambda z: 1 if z == 0 else np.sin(z) / z, x_min = -25, x_max = 25)
        scaledsinc.set_color(YELLOW)
        sinc_title = TexMobject("\\displaystyle \\text{sinc}(x)=\\frac{\\sin x}{x}")
        sinc_title.add_background_rectangle()
        sinc_title.shift(2.5 * UP)
        msinc_title = TexMobject("\\displaystyle \\text{sinc}(\\pi x)=\\frac{\\sin \\pi x}{\\pi x}")
        msinc_title.add_background_rectangle()
        msinc_title.shift(2.5 * UP)
        self.play(Write(sinc), Write(sinc_title))
        self.wait()
        self.play(
            scaledsinc.stretch, 1 / PI, {"dim": 0},
            FadeOut(sinc),
        )
        self.play(Transform(sinc_title, msinc_title))
        self.wait()
        test = ParametricFunction(lambda t: np.array([(2 if t == 0 else 1 + np.sin(PI * t) / PI / t) * np.cos(t), (2 if t == 0 else 1 + np.sin(PI * t) / PI / t) * np.sin(t), 0]), t_min = -10, t_max = 10, color = GREEN)
        test.shift(2.5 * DOWN)
        self.play(
            self.camera_frame.shift, 1 * DOWN,
            Write(test)
        )
        self.wait()
        xval = ValueTracker(-8)
        t = xval.get_value()
        yarrow = Arrow(DOWN + t * RIGHT, DOWN + t * RIGHT + (0 if t == 0 else np.sin(PI * t) / PI / t) * UP + UP)
        garrow = Arrow(ORIGIN, np.array([(2 if t == 0 else 1 + np.sin(PI * t) / PI / t) * np.cos(t), (2 if t == 0 else 1 + np.sin(PI * t) / PI / t) * np.sin(t), 0])).shift(2.5 * DOWN)
        arrow_grp = VGroup(yarrow, garrow)
        self.play(Write(arrow_grp))
        self.wait()
        def arrow_updater(obj):
            y, g = obj
            t = xval.get_value()
            y.become(Arrow(DOWN + t * RIGHT, t * RIGHT + (1 if t == 0 else np.sin(PI * t) / PI / t) * UP))
            g.become(Arrow(ORIGIN, np.array([(2 if t == 0 else 1 + np.sin(PI * t) / PI / t) * np.cos(t), (2 if t == 0 else 1 + np.sin(PI * t) / PI / t) * np.sin(t), 0])).shift(2.5 * DOWN))
        arrow_grp.add_updater(arrow_updater)
        self.add(arrow_grp)
        self.play(xval.set_value, 8, run_time = 25, rate_func = linear)
        self.play(xval.set_value, 0, run_time = 12, rate_func = linear)
        arrow_grp.clear_updaters()
        self.wait(5)

class SmoothingFunction(GraphScene, MovingCameraScene):
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
        self.add(nplane, planegrid)
        self.wait()
        truncexp = self.get_graph(lambda z: np.exp(-abs(z)), x_min = -4, x_max = 4, step_size = 1 / 500)
        fx = TexMobject("f(x)=e^{-|x|}").scale(0.5)
        fx.add_background_rectangle()
        fx.move_to(2.5 * RIGHT + 0.5 * DOWN)
        scaleval = 0.5
        self.play(
            self.camera_frame.scale, scaleval,
            self.camera_frame.shift, 0.5 * UP,
            Write(truncexp),
            Write(fx)
        )
        self.wait()
        title_grp = VGroup(
            #TexMobject("\\displaystyle h_1(x)=\\frac{1}{2}\\left(f\\left(x - \\frac{1}{2}\\right) + f\\left(x + \\frac{1}{2}\\right)\\right)").scale(scaleval).shift(2 * UP).add_background_rectangle(),
            TexMobject("\\displaystyle h_1(x)=\\frac{1}{2}\\left(f\\left(x - 0.5\\right) + f\\left(x + 0.5\\right)\\right)").scale(scaleval).shift(2 * UP).add_background_rectangle(),
            TexMobject("\\displaystyle h_2(x)=\\frac{1}{3}\\left(f\\left(x - 0.5\\right) + f(x) + f\\left(x + 0.5\\right)\\right)").scale(scaleval).shift(2 * UP).add_background_rectangle(),
            TexMobject("\\displaystyle h_N(x)=\\frac{1}{N+1}\\sum_{k=0}^N f\\left(x - \\frac{1}{2} + \\frac{k}{N} \\right)").scale(scaleval).shift(2 * UP).add_background_rectangle(),
            TexMobject("\\displaystyle h(x)=h_\\infty(x)=\\int_{0}^{1}f(x-0.5+t)\\,dt").scale(scaleval).shift(2 * UP).add_background_rectangle(),
        )
        self.play(Write(title_grp[0]))
        self.wait()
        x = 0.25
        xdot = Dot(x * RIGHT, radius = DEFAULT_DOT_RADIUS * scaleval)
        xline = Line(LEFT / 2, RIGHT / 2).shift(x * RIGHT)
        lend, rend = self.input_to_graph_point(x - 0.5, truncexp), self.input_to_graph_point(x + 0.5, truncexp)
        avs = Dot((lend + rend) / 2, radius = DEFAULT_DOT_RADIUS * scaleval)
        avsline = Line(x * RIGHT, (lend + rend) / 2)
        xdivs = VGroup(DashedLine((x - 0.5) * RIGHT, lend), DashedLine((x + 0.5) * RIGHT, rend))
        self.play(Write(xdot), Write(xline))
        self.wait()
        self.play(Write(xdivs))
        self.wait()
        self.play(
            ReplacementTransform(xdivs.copy(), avsline),
            Write(avs),
        )
        self.wait()
        self.play(*[FadeOut(obj) for obj in [xdot, xline, avs, avsline, xdivs]])
        self.wait()
        def dconv_animator(func, numdivs = 2, rtime = 10, xmi = -3.5, xma = 4, colorpref = GREEN):
            x = ValueTracker(xmi)
            def group_creator(y, divs = numdivs):
                ydot = Dot(y * RIGHT, radius = DEFAULT_DOT_RADIUS * scaleval)
                yline = Line(LEFT / 2, RIGHT / 2).shift(y * RIGHT)
                ydivs = VGroup()
                ysum = 0
                for i in range(1 + divs):
                    temp = self.input_to_graph_point(y - 0.5 + i / divs, func)
                    ysum += temp[1]
                    ydivs.add(DashedLine((y - 0.5 + i / divs) * RIGHT, temp))
                ydivs.fade(0.5)
                ysum /= (1 + divs)
                ysum = np.array([y, ysum, 0])
                yavs = Dot(ysum, radius = DEFAULT_DOT_RADIUS * scaleval)
                yavsline = Line(y * RIGHT, ysum)
                return VGroup(ydot, yline, ydivs, yavs, yavsline)
            path = VMobject(color = colorpref)
            path.set_points_as_corners([x.get_value() * RIGHT, x.get_value() * RIGHT + 0.001 * UP])
            conv_grp = VGroup(path, group_creator(x.get_value()))
            def path_updater(obj):
                p, g = obj
                g.become(group_creator(x.get_value()))
                op = p.copy()
                op.append_vectorized_mobject(Line(op.points[-1], g[-2].get_center()))
                p.become(op)
            conv_grp.add_updater(path_updater)
            self.add(conv_grp)
            self.play(x.set_value, xma, run_time = rtime, rate_func = linear)
            conv_grp.clear_updaters()
            self.remove(conv_grp)
            return conv_grp[0]
        for i in range(1, 1 + 3):
            temp = dconv_animator(truncexp, numdivs = i)
            self.add(temp)
            self.play(
                FadeOut(temp),
                FadeInFrom(title_grp[i], UP),
                FadeOutAndShift(title_grp[i - 1], DOWN),
            )
        def convolvedgrph(obj1, obj2):
            #make sure the x-coordinates of both the objects are the same
            #make sure the x-coordinates are symmetric i.e. from +x0 to -x0 for some x0
            ob1 = [x[1] for x in obj1]
            ob2 = [x[1] for x in obj2]
            l = len(ob2)
            res = signal.convolve(ob1, ob2) / sum(ob2)
            l = (l - 1) // 2
            res = res[l:]
            res = res[:-l]
            res = [[x[0], y, 0] for x, y in zip(obj2, res)]
            return res
        def get_idx(arr, x):
            return np.searchsorted(arr, x)
        def get_from_x_coordinate(obj, x):
            temp = obj.points
            indx = get_idx(temp[:, 0], x)
            return temp[indx]
        def area_creator(obj1, x1, x2):
            resobj = VMobject(fill_color = YELLOW, fill_opacity = 0.5)
            tpts = obj1.points.copy()
            idx1, idx2 = get_idx(tpts[:, 0], x1), get_idx(tpts[:, 0], x2)
            tpts = tpts[idx1: idx2]
            tpts = np.insert(tpts, 0, [[x1, 0, 0]], axis = 0)
            tpts = np.append(tpts, [[x2, 0, 0], [x1, 0, 0]], axis = 0)
            resobj.set_points_as_corners(tpts)
            return resobj
        def uconv_animator(funcpts, rtime = 10, colorpref = GREEN):
            eps = funcpts[1][0] - funcpts[0][0]
            xmi, xma = funcpts[0][0], funcpts[-1][0]
            unitbox = [np.array([x, 1 if -0.5 <= x <= 0.5 else 0, 0]) for x in np.arange(xmi, xma + eps, eps)]
            func_g = VMobject().set_points_as_corners(funcpts)
            movsq = VGroup(
            Line(LEFT / 2, LEFT / 2 + UP),
            Line(LEFT / 2 + UP, RIGHT / 2 + UP),
            Line(RIGHT / 2 + UP, RIGHT / 2),
            DashedLine(ORIGIN, UP),
            )
            convpts = convolvedgrph(funcpts, unitbox)
            convpts_g = VMobject().set_points_as_corners(convpts)
            movsq.shift(xmi * RIGHT)
            func_area = area_creator(func_g, movsq.get_critical_point(DOWN)[0] - 0.5, movsq.get_critical_point(DOWN)[0] + 0.5)
            func_path = VMobject(color = colorpref)
            func_dot = Dot(movsq[-1].get_start(), color = colorpref, radius = DEFAULT_DOT_RADIUS / 3)
            func_path.set_points_as_corners([func_dot.get_center(), func_dot.get_center() + 0.001 * UP])
            func_path_grp = VGroup(func_path, func_dot, convpts_g.fade(1))
            func_grp = VGroup(movsq, func_area, func_g.fade(1))
            def conv_updater(obj):
                tm, ta, tc = obj
                curr = tm.get_critical_point(DOWN)[0]
                ta.become(area_creator(tc, curr - 0.5, curr + 0.5))
            def path_updater(obj):
                pat, pd, pc = obj
                curr = get_from_x_coordinate(pc, movsq[-1].get_start()[0])
                pd.move_to(curr)
                old_path = pat.copy()
                old_path.append_vectorized_mobject(Line(old_path.points[-1], curr))
                pat.become(old_path)
            func_grp.add_updater(conv_updater)
            func_path_grp.add_updater(path_updater)
            self.add(func_grp, func_path_grp)
            self.play(movsq.shift, (xma - xmi) * RIGHT, run_time = rtime, rate_func = linear)
            func_grp.clear_updaters()
            func_path_grp.clear_updaters()
            self.remove(func_grp, func_path_grp)
            return convpts
        fpts = [np.array([x, np.exp(-abs(x)), 0]) for x in np.arange(-4, 4 + 1 / 50, 1 / 50)]
        temp = uconv_animator(fpts)
        temp_g = VMobject(color = GREEN).set_points_as_corners(temp)
        self.add(temp_g)
        self.wait()
        self.play(*[FadeOut(obj) for obj in [fx, nplane, planegrid, truncexp, temp_g]])
        self.wait(5)

class ConvolutionDefined(GraphScene, MovingCameraScene):
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
        scaleval = 0.5
        nplane = VGroup(DoubleArrow(8 * LEFT, 8 * RIGHT), DoubleArrow(5.5 * DOWN, 5.5 * UP)).set_color(GREY)
        planegrid = VGroup()
        for i in range(-10, 10 + 1):
            planegrid.add(Line(5 * UP + i * RIGHT, 5 * DOWN + i * RIGHT))
        for i in range(-5, 5 + 1):
            planegrid.add(Line(i * UP + 10 * RIGHT, i * UP + 10 * LEFT))
        planegrid.set_color(GREY)
        prev_eqn = TexMobject("\\displaystyle h(x)=h_\\infty(x)=\\int_{0}^{1}f(x-0.5+t)\\,dt").scale(scaleval).shift(2 * UP).add_background_rectangle()
        planegrid.fade(0.875)
        self.camera_frame.scale(scaleval)
        self.camera_frame.shift(0.5 * UP)
        self.add(prev_eqn)
        self.wait()
        self.play(Transform(prev_eqn, TexMobject("\\displaystyle h(x)=\\int_{0}^{1}f(x-0.5+t)\\,dt").scale(scaleval).shift(2 * UP).add_background_rectangle()))
        self.wait()
        eqn2 = TexMobject("\\displaystyle h(x)=\\int_{x-0.5}^{x+0.5}f(t)\\,dt").scale(scaleval).shift(1 * UP).add_background_rectangle()
        self.play(ReplacementTransform(prev_eqn.copy(), eqn2))
        self.wait()
        eqn3 = TexMobject("\\displaystyle h(x)=\\int\\limits_{-\\infty}^{\\infty}f(t)\\cdot\\text{rect}(x - t)\\,dt").scale(scaleval).shift(0 * UP).set_color(BLUE).add_background_rectangle()
        self.wait()
        rectdefn = TexMobject("\\displaystyle \\text{rect}(x)=\\begin{cases}1, & |x| \\leq 0.5 \\\\ 0, & \\text{otherwise} \\end{cases}").scale(scaleval).shift(-1 * UP).add_background_rectangle()
        self.play(
            FadeInFrom(eqn3, UP),
            FadeInFrom(rectdefn, UP),
        )
        self.wait()
        self.play(
            FadeOut(prev_eqn),
            FadeOut(eqn2),
            FadeOut(rectdefn),
            Transform(eqn3, TexMobject("\\displaystyle h(x)=\\int\\limits_{-\\infty}^{\\infty}f(t)\\cdot g(x - t)\\,dt").scale(scaleval).set_color(BLUE).add_background_rectangle()),
            self.camera_frame.move_to, ORIGIN,
        )
        self.wait()
        scatitle = TextMobject("scaling function").scale(scaleval).move_to(1.5 * RIGHT + 1 * UP).set_color(YELLOW)
        scaarrow = Arrow(0.75 * RIGHT + 0.15 * UP, scatitle.get_critical_point(DOWN), max_tip_length_to_length_ratio = 0.25 * scaleval)
        scagp = VGroup(scatitle, scaarrow)
        self.play(Write(scagp))
        self.wait(5)

class ConvolutionFormula(GraphScene, MovingCameraScene):
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
        scaleval = 0.5
        nplane = VGroup(DoubleArrow(8 * LEFT, 8 * RIGHT), DoubleArrow(5.5 * DOWN, 5.5 * UP)).set_color(GREY)
        planegrid = VGroup()
        for i in range(-10, 10 + 1):
            planegrid.add(Line(5 * UP + i * RIGHT, 5 * DOWN + i * RIGHT))
        for i in range(-5, 5 + 1):
            planegrid.add(Line(i * UP + 10 * RIGHT, i * UP + 10 * LEFT))
        planegrid.set_color(GREY)
        planegrid.fade(0.875)
        self.camera_frame.scale(scaleval)
        self.camera_frame.shift(0.5 * UP)
        rectdefn = TexMobject("\\displaystyle \\text{rect}(t)=\\begin{cases}1, & |t| \\leq 0.5 \\\\ 0, & \\text{otherwise} \\end{cases}").scale(scaleval).shift(2 * UP).add_background_rectangle()
        self.add(nplane, planegrid, rectdefn)
        self.wait()
        rectgrph = VGroup(
            Line(5 * LEFT, LEFT / 2),
            Line(LEFT / 2, LEFT / 2 + UP),
            Line(LEFT / 2 + UP, RIGHT / 2 + UP),
            Line(RIGHT / 2 + UP, RIGHT / 2),
            Line(RIGHT / 2, RIGHT / 2 + 5 * RIGHT),
        )
        rectgrph.set_color(BLUE)
        self.play(ShowCreation(rectgrph))
        self.wait()
        fourect = TexMobject("\\mathcal{F}\\{\\text{rect}(t)\\}").scale(scaleval).add_background_rectangle()
        eqsign = TexMobject("=").scale(scaleval).add_background_rectangle()
        eqsign.shift(0.5 * UP)
        fourect.move_to(0.1 * LEFT + eqsign.get_critical_point(LEFT) + fourect.get_center() - fourect.get_critical_point(RIGHT))
        integ = VGroup(
            TexMobject("\\displaystyle \\int\\limits_{-\\infty}^\\infty \\text{rect}(t)e^{-2\\pi i t \\xi}\\,dt").scale(scaleval).add_background_rectangle(),
            TexMobject("\\displaystyle \\int\\limits_{-1/2}^{1/2} 1\\cdot e^{-2\\pi i t \\xi}\\,dt").scale(scaleval).add_background_rectangle(),
            TexMobject("\\displaystyle \\frac{e^{\\pi i t \\xi}-e^{-\\pi i t \\xi}}{2i\\pi\\xi}").scale(scaleval).add_background_rectangle(),
            TexMobject("\\displaystyle \\frac{\\sin(\\pi\\xi)}{\\pi\\xi}").scale(scaleval).add_background_rectangle()
        )
        for text in integ:
            text.move_to(eqsign.get_critical_point(RIGHT) + text.get_center() - text.get_critical_point(LEFT))
        for text in integ[-2:]:
            text.shift(0.1 * RIGHT)
        self.play(Write(fourect), Write(eqsign), Write(integ[0]))
        self.wait()
        for i in range(1, len(integ)):
            self.play(
                FadeOutAndShift(integ[i - 1], DOWN),
                FadeInFrom(integ[i], UP)
            )
            self.wait()
        sincgrph = self.get_graph(lambda t: 1 if t == 0 else np.sin(PI * t) / PI / t, x_min = -4, x_max = 4).set_color(GREEN)
        eqgrp = VGroup(fourect, eqsign, integ[-1])
        self.play(
            ShowCreation(sincgrph),
            FadeOutAndShift(rectdefn, UP),
            eqgrp.shift, 1.5 * UP,
        )
        self.wait()
        transinteg = VGroup(
            #TextMobject("If "),
            TexMobject("\\mathcal{F}\\{f(t)\\}=\\hat{f}(\\xi)=\\int\\limits_{-\\infty}^{\\infty}e^{-2\\pi i t \\xi}f(t)\\,dt"),
            #TextMobject(","),
            #TextMobject("then "),
            TexMobject("\\mathcal{F}^{-1}\\{\\hat{f}(\\xi)\\}=f(t)=\\int\\limits_{-\\infty}^{\\infty}e^{2\\pi i t \\xi}\\hat{f}(\\xi)\\,d\\xi"),
        )
        transinteg[0][0][2:6].set_color(BLUE)
        transinteg[0][0][8:13].set_color(GREEN)
        transinteg[0][0][19:25].set_color(YELLOW)
        transinteg[0][0][25:29].set_color(BLUE)
        transinteg[1][0][4:9].set_color(GREEN)
        transinteg[1][0][11:15].set_color(BLUE)
        transinteg[1][0][21:26].set_color(YELLOW)
        transinteg[1][0][26:31].set_color(GREEN)
        #transinteg[1].set_color_by_tex("t", GREEN, substring = False)
        #transinteg[0][0].set_color_by_tex("\\hat{f}(\\xi)", GREEN)
        for text in transinteg:
            text.scale_in_place(0.375)
        transinteg[1].move_to(0.125 * DOWN + transinteg[0].get_critical_point(DOWN) + transinteg[1].get_center() - transinteg[1].get_critical_point(UP))
        transinteg.add_background_rectangle()
        transinteg.move_to(0.7 * DOWN)
        self.play(Write(transinteg))
        self.wait()
        invtrans = TexMobject("\\displaystyle \\mathcal{F}\\left\\{\\frac{\\sin(\\pi t)}{\\pi t}\\right\\}=\\text{rect}(\\xi)").scale(scaleval).add_background_rectangle()
        invtrans.move_to(1.75 * RIGHT + 1.5 * UP)
        self.play(
            Write(invtrans),
            eqgrp.move_to, 1.75 * LEFT + 1.5 * UP,
        )
        self.wait()
        atxi = VGroup(
            TexMobject("\\displaystyle \\mathcal{F}\\left\\{\\frac{\\sin(\\pi t)}{\\pi t}\\right\\}_{\\xi = 0}"),
            TexMobject("="),
            TexMobject("\\text{rect}(0)"),
            TexMobject("\\displaystyle \\int\\limits_{-\\infty}^\\infty \\frac{\\sin (\\pi t)}{\\pi t}\\,dt"),
            TexMobject("1"),
        )
        for text in atxi:
            text.scale(scaleval).add_background_rectangle()
        atxi[1].shift(0.5 * DOWN)
        atxi[0].move_to(atxi[1].get_critical_point(LEFT) + atxi[0].get_center() - atxi[0].get_critical_point(RIGHT))
        atxi[2].move_to(0.125 * RIGHT + atxi[1].get_critical_point(RIGHT) + atxi[2].get_center() - atxi[2].get_critical_point(LEFT))
        atxi[3].move_to(0.125 * LEFT + atxi[1].get_critical_point(LEFT) + atxi[3].get_center() - atxi[3].get_critical_point(RIGHT))
        atxi[4].move_to(0.125 * RIGHT + atxi[1].get_critical_point(RIGHT) + atxi[4].get_center() - atxi[4].get_critical_point(LEFT))
        self.play(
            FadeInFrom(atxi[:3], UP),
            FadeOutAndShift(transinteg, DOWN),
        )
        self.wait()
        self.play(
            FadeOutAndShift(atxi[0], DOWN),
            FadeInFrom(atxi[3], UP),
            FadeOutAndShift(atxi[2], UP),
            FadeInFrom(atxi[4], DOWN),
        )
        self.play(
            VGroup(atxi[1], atxi[3], atxi[4]).move_to, 1.75 * UP,
            eqgrp.shift, 2.25 * DOWN,
            invtrans.shift, 2.25 * DOWN,
        )            
        self.wait(5)

class ConvolutionTheorem(GraphScene, MovingCameraScene):
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
        scaleval = 0.5
        nplane = VGroup(DoubleArrow(8 * LEFT, 8 * RIGHT), DoubleArrow(5.5 * DOWN, 5.5 * UP)).set_color(GREY)
        planegrid = VGroup()
        for i in range(-10, 10 + 1):
            planegrid.add(Line(5 * UP + i * RIGHT, 5 * DOWN + i * RIGHT))
        for i in range(-5, 5 + 1):
            planegrid.add(Line(i * UP + 10 * RIGHT, i * UP + 10 * LEFT))
        planegrid.set_color(GREY)
        planegrid.fade(0.875)
        self.camera_frame.scale(scaleval)
        titl = TexMobject("\\underline{\\text{Convolution Theorem}}").scale(0.625).set_color(YELLOW).shift(1.75 * UP)
        self.play(Write(titl))
        self.wait()
        convtxt = VGroup(
            TexMobject("\\mathcal{F}\\{f*g\\}=\\mathcal{F}\\{f\\} \\cdot \\mathcal{F}\\{g\\}"),
            TexMobject("\\mathcal{F}\\{f \\cdot g\\}=\\mathcal{F}\\{f\\} * \\mathcal{F}\\{g\\}"),
        )
        for txt in convtxt:
            txt.scale_in_place(scaleval)
        convtxt[0].shift(0.75 * UP)
        convtxt[1].shift(0.75 * DOWN)
        self.play(Write(convtxt))
        self.wait(5)
        
class SincSquared(GraphScene, MovingCameraScene):
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
        scaleval = 0.5
        nplane = VGroup(DoubleArrow(8 * LEFT, 8 * RIGHT), DoubleArrow(5.5 * DOWN, 5.5 * UP)).set_color(GREY)
        planegrid = VGroup()
        for i in range(-10, 10 + 1):
            planegrid.add(Line(5 * UP + i * RIGHT, 5 * DOWN + i * RIGHT))
        for i in range(-5, 5 + 1):
            planegrid.add(Line(i * UP + 10 * RIGHT, i * UP + 10 * LEFT))
        planegrid.set_color(GREY)
        planegrid.fade(0.875)
        self.camera_frame.scale(scaleval)
        self.camera_frame.shift(0.5 * UP)
        sincsq = TexMobject("\\displaystyle \\int\\limits_{-\\infty}^\\infty \\text{sinc}^2(\\pi t)\\,dt").scale(scaleval).add_background_rectangle()
        sincsq.shift(1.5 * UP)
        self.play(Write(sincsq))
        self.wait()
        sincsqf = VGroup(
            TexMobject("\\mathcal{F}\\{\\text{sinc}^2(\\pi t)\\}"),
            TexMobject("=\\displaystyle \\int\\limits_{-\\infty}^\\infty e^{-2\\pi i t \\xi}\\text{sinc}^2(\\pi t)\\,dt"),
            TexMobject("=\\mathcal{F}\\{\\text{sinc}(\\pi t) \\cdot \\text{sinc}(\\pi t)\\}"),
            TexMobject("=\\mathcal{F}\\{\\text{sinc}(\\pi t)\\} * \\mathcal{F}\\{\\text{sinc}(\\pi t)\\}"),
            TexMobject("=\\text{rect}(\\xi) * \\text{rect}(\\xi)"),
        )
        for text in sincsqf:
            text.scale_in_place(scaleval)
            text.add_background_rectangle()
        sincsqf[0].shift(1 * LEFT)
        for text in sincsqf[1:]:
            text.next_to(sincsqf[0], RIGHT, DEFAULT_MOBJECT_TO_MOBJECT_BUFFER * scaleval)
        self.play(FadeIn(sincsqf[0]), FadeIn(sincsqf[1]))
        self.wait()
        for i in range(2, len(sincsqf)):
            self.play(
                FadeOutAndShift(sincsqf[i - 1], DOWN),
                FadeInFrom(sincsqf[i], UP)
            )
            self.wait()
        eqngp = VGroup(sincsqf[0], sincsqf[-1])
        self.play(
            FadeOut(sincsq),
            eqngp.shift, 1.5 * UP,
            FadeIn(nplane),
            FadeIn(planegrid),
        )
        self.add(eqngp)
        self.wait()
        boxgraph = VGroup(
            Line(5 * LEFT, LEFT / 2),
            Line(LEFT / 2, LEFT / 2 + UP),
            Line(LEFT / 2 + UP, RIGHT / 2 + UP),
            Line(RIGHT / 2 + UP, RIGHT / 2),
            Line(RIGHT / 2, RIGHT / 2 + 5 * RIGHT),
        )
        boxgraph.set_color(BLUE)
        self.play(Write(boxgraph))
        self.wait()
        def convolvedgrph(obj1, obj2):
            #make sure the x-coordinates of both the objects are the same
            ob1 = [x[1] for x in obj1]
            ob2 = [x[1] for x in obj2]
            l = len(ob2)
            res = signal.convolve(ob1, ob2) / sum(ob2)
            l = (l - 1) // 2
            res = res[l:]
            res = res[:-l]
            res = [[x[0], y, 0] for x, y in zip(obj2, res)]
            return res
        def get_idx(arr, x):
            return np.searchsorted(arr, x)
        def get_from_x_coordinate(obj, x):
            temp = obj.points
            indx = get_idx(temp[:, 0], x)
            return temp[indx]
        def area_creator(obj1, x1, x2):
            resobj = VMobject(fill_color = YELLOW, fill_opacity = 0.5)
            tpts = obj1.points.copy()
            idx1, idx2 = get_idx(tpts[:, 0], x1), get_idx(tpts[:, 0], x2)
            tpts = tpts[idx1: idx2]
            tpts = np.insert(tpts, 0, [[x1, 0, 0]], axis = 0)
            tpts = np.append(tpts, [[x2, 0, 0], [x1, 0, 0]], axis = 0)
            resobj.set_points_as_corners(tpts)
            return resobj
        def uconv_animator(funcpts, rtime = 10, colorpref = GREEN):
            eps = funcpts[1][0] - funcpts[0][0]
            xmi, xma = funcpts[0][0], funcpts[-1][0]
            unitbox = [np.array([x, 1 if -0.5 <= x <= 0.5 else 0, 0]) for x in np.arange(xmi, xma + eps, eps)]
            func_g = VMobject().set_points_as_corners(funcpts)
            movsq = VGroup(
            Line(LEFT / 2, LEFT / 2 + UP),
            Line(LEFT / 2 + UP, RIGHT / 2 + UP),
            Line(RIGHT / 2 + UP, RIGHT / 2),
            DashedLine(ORIGIN, UP),
            )
            convpts = convolvedgrph(funcpts, unitbox)
            convpts_g = VMobject().set_points_as_corners(convpts)
            movsq.shift(xmi * RIGHT)
            func_area = area_creator(func_g, movsq.get_critical_point(DOWN)[0] - 0.5, movsq.get_critical_point(DOWN)[0] + 0.5)
            func_path = VMobject(color = colorpref)
            func_dot = Dot(movsq[-1].get_start(), color = colorpref, radius = DEFAULT_DOT_RADIUS / 3)
            func_path.set_points_as_corners([func_dot.get_center(), func_dot.get_center() + 0.001 * UP])
            func_path_grp = VGroup(func_path, func_dot, convpts_g.fade(1))
            func_grp = VGroup(movsq, func_area, func_g.fade(1))
            def conv_updater(obj):
                tm, ta, tc = obj
                curr = tm.get_critical_point(DOWN)[0]
                ta.become(area_creator(tc, curr - 0.5, curr + 0.5))
            def path_updater(obj):
                pat, pd, pc = obj
                curr = get_from_x_coordinate(pc, movsq[-1].get_start()[0])
                pd.move_to(curr)
                old_path = pat.copy()
                old_path.append_vectorized_mobject(Line(old_path.points[-1], curr))
                pat.become(old_path)
            func_grp.add_updater(conv_updater)
            func_path_grp.add_updater(path_updater)
            self.add(func_grp, func_path_grp)
            self.play(movsq.shift, (xma - xmi) * RIGHT, run_time = rtime, rate_func = linear)
            func_grp.clear_updaters()
            func_path_grp.clear_updaters()
            self.remove(func_grp, func_path_grp)
            return convpts
        fpts = [np.array([x, 1 if -0.5 <= x <= 0.5 else 0, 0]) for x in np.arange(-4, 4 + 1 / 75, 1 / 75)]
        temp = uconv_animator(fpts)
        temp_g = VMobject(color = GREEN).set_points_as_corners(temp)
        self.add(temp_g)
        self.wait()
        tritext = TexMobject("=\\text{tri}(\\xi)").scale(scaleval)
        tritext.next_to(sincsqf[0], RIGHT, DEFAULT_MOBJECT_TO_MOBJECT_BUFFER * scaleval)
        self.play(FadeOutAndShift(sincsqf[-1], DOWN), FadeInFrom(tritext, UP))
        self.wait()
        iarea = area_creator(temp_g, -0.5, 0.5)
        feqn = VGroup(
            TexMobject("\\mathcal{F}\\{\\text{sinc}^2(\\pi t)\\}_{\\xi=0}"),
            TexMobject("="),
            TexMobject("\\text{tri}(0)"),
        )
        for text in feqn:
            text.scale_in_place(scaleval)
            text.add_background_rectangle()
        feqn.shift(1.5 * UP + 0.5 * RIGHT)
        feqn[0].next_to(feqn[1], LEFT, DEFAULT_MOBJECT_TO_MOBJECT_BUFFER * scaleval)
        feqn[2].next_to(feqn[1], RIGHT, DEFAULT_MOBJECT_TO_MOBJECT_BUFFER * scaleval)
        self.play(FadeOutAndShift(sincsqf[0], RIGHT), FadeOutAndShift(tritext, RIGHT), FadeInFrom(feqn, LEFT), Write(iarea))
        self.wait()
        feqnu = VGroup(
            TexMobject("\\int\\limits_{-\\infty}^\\infty \\text{sinc}^2(\\pi t)\\,dt"),
            TexMobject("="),
            TexMobject("1"),
        )
        for text in feqnu:
            text.scale_in_place(scaleval)
            text.add_background_rectangle()
        feqnu.shift(1.5 * UP + 0.5 * RIGHT)
        feqnu[0].next_to(feqnu[1], LEFT, DEFAULT_MOBJECT_TO_MOBJECT_BUFFER * scaleval)
        feqnu[2].next_to(feqnu[1], RIGHT, DEFAULT_MOBJECT_TO_MOBJECT_BUFFER * scaleval)
        self.play(
            FadeOutAndShift(feqn[0], DOWN), FadeInFrom(feqnu[0], UP),
            FadeOutAndShift(feqn[2], UP), FadeInFrom(feqnu[2], DOWN),
        )
        self.wait()
        self.play(
            FadeOut(iarea),
            FadeOut(boxgraph),
            temp_g.set_color, BLUE,
            FadeOut(feqnu[0]), FadeOut(feqnu[2]), FadeOut(feqn[1]),
        )
        self.wait(5)

class SincCubed(GraphScene, MovingCameraScene):
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
        scaleval = 0.5
        nplane = VGroup(DoubleArrow(8 * LEFT, 8 * RIGHT), DoubleArrow(5.5 * DOWN, 5.5 * UP)).set_color(GREY)
        planegrid = VGroup()
        for i in range(-10, 10 + 1):
            planegrid.add(Line(5 * UP + i * RIGHT, 5 * DOWN + i * RIGHT))
        for i in range(-5, 5 + 1):
            planegrid.add(Line(i * UP + 10 * RIGHT, i * UP + 10 * LEFT))
        planegrid.set_color(GREY)
        planegrid.fade(0.875)
        self.camera_frame.scale(scaleval)
        self.camera_frame.shift(0.5 * UP)
        fpts = [np.array([x, 1 if -0.5 <= x <= 0.5 else 0, 0]) for x in np.arange(-4, 4 + 1 / 75, 1 / 75)]
        def convolvedgrph(obj1, obj2):
            #make sure the x-coordinates of both the objects are the same
            ob1 = [x[1] for x in obj1]
            ob2 = [x[1] for x in obj2]
            l = len(ob2)
            res = signal.convolve(ob1, ob2) / sum(ob2)
            l = (l - 1) // 2
            res = res[l:]
            res = res[:-l]
            res = [[x[0], y, 0] for x, y in zip(obj2, res)]
            return res
        def get_idx(arr, x):
            return np.searchsorted(arr, x)
        def get_from_x_coordinate(obj, x):
            temp = obj.points
            indx = get_idx(temp[:, 0], x)
            return temp[indx]
        def area_creator(obj1, x1, x2):
            resobj = VMobject(fill_color = YELLOW, fill_opacity = 0.5)
            tpts = obj1.points.copy()
            idx1, idx2 = get_idx(tpts[:, 0], x1), get_idx(tpts[:, 0], x2)
            tpts = tpts[idx1: idx2]
            tpts = np.insert(tpts, 0, [[x1, 0, 0]], axis = 0)
            tpts = np.append(tpts, [[x2, 0, 0], [x1, 0, 0]], axis = 0)
            resobj.set_points_as_corners(tpts)
            return resobj
        def uconv_animator(funcpts, rtime = 10, colorpref = GREEN):
            eps = funcpts[1][0] - funcpts[0][0]
            xmi, xma = funcpts[0][0], funcpts[-1][0]
            unitbox = [np.array([x, 1 if -0.5 <= x <= 0.5 else 0, 0]) for x in np.arange(xmi, xma + eps, eps)]
            func_g = VMobject().set_points_as_corners(funcpts)
            movsq = VGroup(
            Line(LEFT / 2, LEFT / 2 + UP),
            Line(LEFT / 2 + UP, RIGHT / 2 + UP),
            Line(RIGHT / 2 + UP, RIGHT / 2),
            DashedLine(ORIGIN, UP),
            )
            convpts = convolvedgrph(funcpts, unitbox)
            convpts_g = VMobject().set_points_as_corners(convpts)
            movsq.shift(xmi * RIGHT)
            func_area = area_creator(func_g, movsq.get_critical_point(DOWN)[0] - 0.5, movsq.get_critical_point(DOWN)[0] + 0.5)
            func_path = VMobject(color = colorpref)
            func_dot = Dot(movsq[-1].get_start(), color = colorpref, radius = DEFAULT_DOT_RADIUS / 3)
            func_path.set_points_as_corners([func_dot.get_center(), func_dot.get_center() + 0.001 * UP])
            func_path_grp = VGroup(func_path, func_dot, convpts_g.fade(1))
            func_grp = VGroup(movsq, func_area, func_g.fade(1))
            def conv_updater(obj):
                tm, ta, tc = obj
                curr = tm.get_critical_point(DOWN)[0]
                ta.become(area_creator(tc, curr - 0.5, curr + 0.5))
            def path_updater(obj):
                pat, pd, pc = obj
                curr = get_from_x_coordinate(pc, movsq[-1].get_start()[0])
                pd.move_to(curr)
                old_path = pat.copy()
                old_path.append_vectorized_mobject(Line(old_path.points[-1], curr))
                pat.become(old_path)
            func_grp.add_updater(conv_updater)
            func_path_grp.add_updater(path_updater)
            self.add(func_grp, func_path_grp)
            self.play(movsq.shift, (xma - xmi) * RIGHT, run_time = rtime, rate_func = linear)
            func_grp.clear_updaters()
            func_path_grp.clear_updaters()
            self.remove(func_grp, func_path_grp)
            return convpts
        fpts = [np.array([x, 1 if -0.5 <= x <= 0.5 else 0, 0]) for x in np.arange(-4, 4 + 1 / 75, 1 / 75)]
        tripts = convolvedgrph(fpts, fpts)
        trigrph = VMobject(color = BLUE).set_points_as_corners(tripts)
        self.add(nplane, planegrid, trigrph)
        self.wait()
        cubetxt = VGroup(
            TexMobject("\\mathcal{F}\\{\\text{sinc}^3(\\pi t)\\}"),
            TexMobject("="),
            TexMobject("\\mathcal{F}\\{\\text{sinc}^2(\\pi t) \\cdot \\text{sinc}(\\pi t) \\}"),
            TexMobject("\\mathcal{F}\\{\\text{sinc}^2(\\pi t)\\}*\\mathcal{F}\\{\\text{sinc}(\\pi t)\\}"),
            TexMobject("\\text{tri}(\\xi)*\\text{rect}(\\xi)"),
        )
        for text in cubetxt:
            text.scale_in_place(scaleval)
            text.add_background_rectangle()
        cubetxt[1].move_to(1.5 * UP)
        cubetxt[0].next_to(cubetxt[1], LEFT, DEFAULT_MOBJECT_TO_MOBJECT_BUFFER * scaleval)
        for txt in cubetxt[2:]:
            txt.next_to(cubetxt[1], RIGHT, DEFAULT_MOBJECT_TO_MOBJECT_BUFFER * scaleval)
        self.play(Write(cubetxt[:3]))
        self.wait()
        for i in range(3, len(cubetxt)):
            self.play(FadeOutAndShift(cubetxt[i - 1], DOWN), FadeInFrom(cubetxt[i], UP))
            self.wait()
        temp = uconv_animator(tripts)
        temp_g = VMobject(color = GREEN).set_points_as_corners(temp)
        self.add(temp_g)
        self.wait()
        iarea = area_creator(trigrph, -0.5, 0.5)
        feqn = VGroup(
            TexMobject("\\mathcal{F}\\{\\text{sinc}^3(\\pi t)\\}_{\\xi=0}"),
            TexMobject("\\text{tri}(\\xi)*\\text{rect}(\\xi)\\bigg\\rvert_{\\xi = 0}"),
        )
        for text in feqn:
            text.scale_in_place(scaleval)
            text.add_background_rectangle()
        feqn[0].next_to(cubetxt[1], LEFT, DEFAULT_MOBJECT_TO_MOBJECT_BUFFER * scaleval)
        feqn[1].next_to(cubetxt[1], RIGHT, DEFAULT_MOBJECT_TO_MOBJECT_BUFFER * scaleval)
        leqn = TexMobject("\\displaystyle\\int\\limits_{-\\infty}^{\\infty}\\text{sinc}^3(\\pi t)\\,dt=\\frac{3}{4}")
        leqn.scale_in_place(scaleval)
        leqn.add_background_rectangle()
        leqn.shift(1.5 * UP)
        self.play(
            FadeOutAndShift(cubetxt[0], UP), FadeInFrom(feqn[0], DOWN),
            FadeOutAndShift(cubetxt[-1], DOWN), FadeInFrom(feqn[1], UP),
        )
        self.wait()
        self.play(Write(iarea))
        self.wait()
        feqngp = VGroup(feqn[0], cubetxt[1], feqn[1])
        self.play(
            FadeOutAndShift(feqngp, RIGHT), FadeInFrom(leqn, RIGHT),
        )
        self.wait()
        self.play(
            FadeOut(leqn),
            temp_g.set_color, BLUE,
            FadeOut(iarea),
            FadeOut(trigrph)
        )
        self.wait(5)

class SincFourth(GraphScene, MovingCameraScene):
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
        scaleval = 0.5
        nplane = VGroup(DoubleArrow(8 * LEFT, 8 * RIGHT), DoubleArrow(5.5 * DOWN, 5.5 * UP)).set_color(GREY)
        planegrid = VGroup()
        for i in range(-10, 10 + 1):
            planegrid.add(Line(5 * UP + i * RIGHT, 5 * DOWN + i * RIGHT))
        for i in range(-5, 5 + 1):
            planegrid.add(Line(i * UP + 10 * RIGHT, i * UP + 10 * LEFT))
        planegrid.set_color(GREY)
        planegrid.fade(0.875)
        self.camera_frame.scale(scaleval)
        self.camera_frame.shift(0.5 * UP)
        fpts = [np.array([x, 1 if -0.5 <= x <= 0.5 else 0, 0]) for x in np.arange(-4, 4 + 1 / 75, 1 / 75)]
        def convolvedgrph(obj1, obj2):
            #make sure the x-coordinates of both the objects are the same
            ob1 = [x[1] for x in obj1]
            ob2 = [x[1] for x in obj2]
            l = len(ob2)
            res = signal.convolve(ob1, ob2) / sum(ob2)
            l = (l - 1) // 2
            res = res[l:]
            res = res[:-l]
            res = [[x[0], y, 0] for x, y in zip(obj2, res)]
            return res
        def get_idx(arr, x):
            return np.searchsorted(arr, x)
        def get_from_x_coordinate(obj, x):
            temp = obj.points
            indx = get_idx(temp[:, 0], x)
            return temp[indx]
        def area_creator(obj1, x1, x2):
            resobj = VMobject(fill_color = YELLOW, fill_opacity = 0.5)
            tpts = obj1.points.copy()
            idx1, idx2 = get_idx(tpts[:, 0], x1), get_idx(tpts[:, 0], x2)
            tpts = tpts[idx1: idx2]
            tpts = np.insert(tpts, 0, [[x1, 0, 0]], axis = 0)
            tpts = np.append(tpts, [[x2, 0, 0], [x1, 0, 0]], axis = 0)
            resobj.set_points_as_corners(tpts)
            return resobj
        def uconv_animator(funcpts, rtime = 10, colorpref = GREEN):
            eps = funcpts[1][0] - funcpts[0][0]
            xmi, xma = funcpts[0][0], funcpts[-1][0]
            unitbox = [np.array([x, 1 if -0.5 <= x <= 0.5 else 0, 0]) for x in np.arange(xmi, xma + eps, eps)]
            func_g = VMobject().set_points_as_corners(funcpts)
            movsq = VGroup(
            Line(LEFT / 2, LEFT / 2 + UP),
            Line(LEFT / 2 + UP, RIGHT / 2 + UP),
            Line(RIGHT / 2 + UP, RIGHT / 2),
            DashedLine(ORIGIN, UP),
            )
            convpts = convolvedgrph(funcpts, unitbox)
            convpts_g = VMobject().set_points_as_corners(convpts)
            movsq.shift(xmi * RIGHT)
            func_area = area_creator(func_g, movsq.get_critical_point(DOWN)[0] - 0.5, movsq.get_critical_point(DOWN)[0] + 0.5)
            func_path = VMobject(color = colorpref)
            func_dot = Dot(movsq[-1].get_start(), color = colorpref, radius = DEFAULT_DOT_RADIUS / 3)
            func_path.set_points_as_corners([func_dot.get_center(), func_dot.get_center() + 0.001 * UP])
            func_path_grp = VGroup(func_path, func_dot, convpts_g.fade(1))
            func_grp = VGroup(movsq, func_area, func_g.fade(1))
            def conv_updater(obj):
                tm, ta, tc = obj
                curr = tm.get_critical_point(DOWN)[0]
                ta.become(area_creator(tc, curr - 0.5, curr + 0.5))
            def path_updater(obj):
                pat, pd, pc = obj
                curr = get_from_x_coordinate(pc, movsq[-1].get_start()[0])
                pd.move_to(curr)
                old_path = pat.copy()
                old_path.append_vectorized_mobject(Line(old_path.points[-1], curr))
                pat.become(old_path)
            func_grp.add_updater(conv_updater)
            func_path_grp.add_updater(path_updater)
            self.add(func_grp, func_path_grp)
            self.play(movsq.shift, (xma - xmi) * RIGHT, run_time = rtime, rate_func = linear)
            func_grp.clear_updaters()
            func_path_grp.clear_updaters()
            self.remove(func_grp, func_path_grp)
            return convpts
        fpts = [np.array([x, 1 if -0.5 <= x <= 0.5 else 0, 0]) for x in np.arange(-4, 4 + 1 / 75, 1 / 75)]
        tripts = convolvedgrph(fpts, fpts)
        cubepts = convolvedgrph(tripts, fpts)
        trigrph = VMobject(color = BLUE).set_points_as_corners(tripts)
        cubegrph = VMobject(color = BLUE).set_points_as_corners(cubepts)
        self.add(nplane, planegrid, cubegrph)
        self.wait()
        cubetxt = VGroup(
            TexMobject("\\mathcal{F}\\{\\text{sinc}^4(\\pi t)\\}"),
            TexMobject("="),
            TexMobject("\\mathcal{F}\\{\\text{sinc}^3(\\pi t) \\cdot \\text{sinc}(\\pi t) \\}"),
            TexMobject("\\mathcal{F}\\{\\text{sinc}^3(\\pi t)\\}*\\mathcal{F}\\{\\text{sinc}(\\pi t)\\}"),
            TexMobject("\\text{rect}^{*3}(\\xi)*\\text{rect}(\\xi)"),
        )
        for text in cubetxt:
            text.scale_in_place(scaleval)
            text.add_background_rectangle()
        cubetxt[1].move_to(1.5 * UP)
        cubetxt[0].next_to(cubetxt[1], LEFT, DEFAULT_MOBJECT_TO_MOBJECT_BUFFER * scaleval)
        for txt in cubetxt[2:]:
            txt.next_to(cubetxt[1], RIGHT, DEFAULT_MOBJECT_TO_MOBJECT_BUFFER * scaleval)
        self.play(Write(cubetxt[:3]))
        self.wait()
        for i in range(3, len(cubetxt)):
            self.play(FadeOutAndShift(cubetxt[i - 1], DOWN), FadeInFrom(cubetxt[i], UP))
            self.wait()
        temp = uconv_animator(cubepts)
        temp_g = VMobject(color = GREEN).set_points_as_corners(temp)
        self.add(temp_g)
        self.wait()
        iarea = area_creator(cubegrph, -0.5, 0.5)
        feqn = VGroup(
            TexMobject("\\mathcal{F}\\{\\text{sinc}^4(\\pi t)\\}_{\\xi=0}"),
            TexMobject("\\text{rect}^{*3}(\\xi)*\\text{rect}(\\xi)\\bigg\\rvert_{\\xi = 0}"),
        )
        for text in feqn:
            text.scale_in_place(scaleval)
            text.add_background_rectangle()
        feqn[0].next_to(cubetxt[1], LEFT, DEFAULT_MOBJECT_TO_MOBJECT_BUFFER * scaleval)
        feqn[1].next_to(cubetxt[1], RIGHT, DEFAULT_MOBJECT_TO_MOBJECT_BUFFER * scaleval)
        self.play(
            FadeOutAndShift(cubetxt[0], DOWN), FadeInFrom(feqn[0], UP),
            FadeOutAndShift(cubetxt[-1], DOWN), FadeInFrom(feqn[1], UP),
        )
        self.wait()
        self.play(Write(iarea))
        self.wait()
        self.play(
            FadeOutAndShift(feqn[0], UP), FadeInFrom(cubetxt[0], DOWN),
            FadeOutAndShift(feqn[1], UP), FadeInFrom(cubetxt[-1], DOWN),
            FadeOut(iarea),
        )
        self.wait()
        cubetxtm = VGroup(
            TexMobject("\\mathcal{F}\\{\\text{sinc}^4(\\pi t)\\}"),
            TexMobject("="),
            TexMobject("\\mathcal{F}\\{\\text{sinc}^2(\\pi t) \\cdot \\text{sinc}^2(\\pi t) \\}"),
            TexMobject("\\mathcal{F}\\{\\text{sinc}^2(\\pi t)\\}*\\mathcal{F}\\{\\text{sinc}^2(\\pi t)\\}"),
            TexMobject("\\text{tri}(\\xi)*\\text{tri}(\\xi)"),
        )
        for text in cubetxtm:
            text.scale_in_place(scaleval)
            text.add_background_rectangle()
        cubetxtm[1].move_to(1.5 * UP)
        cubetxtm[0].next_to(cubetxtm[1], LEFT, DEFAULT_MOBJECT_TO_MOBJECT_BUFFER * scaleval)
        for txt in cubetxtm[2:]:
            txt.next_to(cubetxtm[1], RIGHT, DEFAULT_MOBJECT_TO_MOBJECT_BUFFER * scaleval)
        self.play(FadeOutAndShift(cubetxt[-1], RIGHT), FadeInFrom(cubetxtm[2], UP))
        self.wait()
        for i in range(3, len(cubetxtm)):
            self.play(FadeOutAndShift(cubetxtm[i - 1], DOWN), FadeInFrom(cubetxtm[i], UP))
            self.wait()
        self.play(FadeOut(temp_g), FadeOut(cubegrph), FadeIn(trigrph))
        self.wait()
        scal = TextMobject("scaling function").scale(scaleval).set_color(YELLOW)
        scal.add_background_rectangle()
        scal.move_to(2.5 * RIGHT + 1 * UP)
        scalarrow = Arrow(1.375 * UP + 1.375 * RIGHT, scal.get_critical_point(UP), max_tip_length_to_length_ratio = 0.25 * scaleval)
        scalgrp = VGroup(scal, scalarrow)
        self.play(Write(scalgrp))
        self.wait()
        textval = 0.375
        rdot = Dot(0.5 * RIGHT + 0.5 * UP, radius = DEFAULT_DOT_RADIUS * textval).set_color(PURPLE)
        rtxt = TexMobject("(t,1-t)").scale(textval).set_color(GREEN).next_to(rdot, RIGHT, 0.1).add_background_rectangle()
        rrdot = Dot(0.5 * RIGHT + 0.25 * UP, radius = DEFAULT_DOT_RADIUS * textval).set_color(PURPLE)
        rrtxt = TexMobject("(t,(1-t)^2)").scale(textval).set_color(GREEN).next_to(rrdot, RIGHT, 0.1).add_background_rectangle()
        self.play(Write(rdot), Write(rtxt))
        self.wait()
        self.play(
            ReplacementTransform(rdot.copy(), rrdot),
            ReplacementTransform(rtxt.copy(), rrtxt),
        )
        self.wait()
        paraarcs = VMobject(fill_color = YELLOW, fill_opacity = 0.5, color = GREEN)
        paraarcspts = [np.array([t, (1 - abs(t)) * (1 - abs(t)), 0]) for t in np.arange(-1, 1 + 1 / 50, 1 / 50)]
        paraarcspts.append(np.array([-1, 0, 0]))
        paraarcs.set_points_as_corners(paraarcspts)
        self.play(
            *[FadeOut(obj) for obj in [rdot, rtxt, rrdot, rrtxt, scalgrp]],
            FadeIn(paraarcs),
        )
        self.wait()
        ffeqn = VGroup(
            TexMobject("\\mathcal{F}\\{\\text{sinc}^4(\\pi t)\\}_{\\xi=0}"),
            TexMobject("\\text{tri}(\\xi)*\\text{tri}(\\xi)\\bigg\\rvert_{\\xi = 0}"),
        )
        for text in ffeqn:
            text.scale_in_place(scaleval)
            text.add_background_rectangle()
        ffeqn[0].next_to(cubetxt[1], LEFT, DEFAULT_MOBJECT_TO_MOBJECT_BUFFER * scaleval)
        ffeqn[1].next_to(cubetxt[1], RIGHT, DEFAULT_MOBJECT_TO_MOBJECT_BUFFER * scaleval)
        leqn = TexMobject("\\displaystyle\\int\\limits_{-\\infty}^{\\infty}\\text{sinc}^4(\\pi t)\\,dt=\\frac{2}{3}")
        leqn.scale_in_place(scaleval)
        leqn.add_background_rectangle()
        leqn.shift(1.5 * UP)
        self.play(
            FadeOutAndShift(cubetxt[0], UP), FadeInFrom(ffeqn[0], DOWN),
            FadeOutAndShift(cubetxtm[-1], DOWN), FadeInFrom(ffeqn[1], UP),
        )
        self.wait()
        ffeqngp = VGroup(ffeqn[0], ffeqn[1], cubetxt[1])
        self.play(
            FadeOutAndShift(ffeqngp, RIGHT), FadeInFrom(leqn, LEFT),
        )
        self.wait(5)

class UserChallenge(GraphScene, MovingCameraScene):
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
        scaleval = 0.5
        nplane = VGroup(DoubleArrow(8 * LEFT, 8 * RIGHT), DoubleArrow(5.5 * DOWN, 5.5 * UP)).set_color(GREY)
        planegrid = VGroup()
        for i in range(-10, 10 + 1):
            planegrid.add(Line(5 * UP + i * RIGHT, 5 * DOWN + i * RIGHT))
        for i in range(-5, 5 + 1):
            planegrid.add(Line(i * UP + 10 * RIGHT, i * UP + 10 * LEFT))
        planegrid.set_color(GREY)
        planegrid.fade(0.875)
        self.camera_frame.scale(scaleval)
        #self.camera_frame.shift(0.5 * UP)
        txtgp = VGroup(
            TexMobject("\\text{For }0 < a \\leq b,"),
            TexMobject("\\displaystyle \\int\\limits_{-\\infty}^{\\infty}\\text{sinc}\\left(\\frac{\\pi t}{a}\\right)\\text{sinc}\\left(\\frac{\\pi t}{b}\\right)\\,dt=a"),
            TexMobject("\\displaystyle \\int\\limits_{-\\infty}^{\\infty}\\text{sinc}^2\\left(\\frac{\\pi t}{a}\\right)\\text{sinc}^2\\left(\\frac{\\pi t}{b}\\right)\\,dt=\\frac{2}{3}\\cdot \\frac{1}{b} \\cdot a^2+a^2\\left(\\frac{1}{a}-\\frac{1}{b}\\right)"),
        )
        for txt in txtgp:
            txt.scale_in_place(scaleval)
        txtgp[0].shift(1.5 * UP)
        txtgp[1].shift(0.5 * UP)
        txtgp[2].shift(0.5 * DOWN)
        for txt in txtgp:
            self.play(Write(txt))
            self.wait()
        self.play(FadeOut(txtgp))
        self.wait()
        borintegs = VGroup(
            TexMobject("\\displaystyle \\int\\limits_{-\\infty}^{\\infty}\\text{sinc}(\\pi t)\\,dt=1"),
            TexMobject("\\displaystyle \\int\\limits_{-\\infty}^{\\infty}\\text{sinc}(\\pi t)\\text{sinc}\\left(\\frac{\\pi t}{3}\\right)\\,dt=1"),
            TexMobject("\\displaystyle \\int\\limits_{-\\infty}^{\\infty}\\text{sinc}(\\pi t)\\text{sinc}\\left(\\frac{\\pi t}{3}\\right)\\text{sinc}\\left(\\frac{\\pi t}{5}\\right)\\,dt=1"),
            TexMobject("\\displaystyle \\int\\limits_{-\\infty}^{\\infty}\\text{sinc}(\\pi t)\\text{sinc}\\left(\\frac{\\pi t}{3}\\right)\\text{sinc}\\left(\\frac{\\pi t}{5}\\right)\\text{sinc}\\left(\\frac{\\pi t}{7}\\right)\\,dt=1"),
            TexMobject("\\displaystyle \\int\\limits_{-\\infty}^{\\infty}\\text{sinc}(\\pi t)\\text{sinc}\\left(\\frac{\\pi t}{3}\\right)\\text{sinc}\\left(\\frac{\\pi t}{5}\\right)\\text{sinc}\\left(\\frac{\\pi t}{7}\\right)\\text{sinc}\\left(\\frac{\\pi t}{9}\\right)\\,dt=1"),
            TexMobject("\\displaystyle \\int\\limits_{-\\infty}^{\\infty}\\text{sinc}(\\pi t)\\text{sinc}\\left(\\frac{\\pi t}{3}\\right)\\text{sinc}\\left(\\frac{\\pi t}{5}\\right)\\text{sinc}\\left(\\frac{\\pi t}{7}\\right)\\text{sinc}\\left(\\frac{\\pi t}{9}\\right)\\text{sinc}\\left(\\frac{\\pi t}{11}\\right)\\,dt=1"),
            TexMobject("\\displaystyle \\int\\limits_{-\\infty}^{\\infty}\\text{sinc}(\\pi t)\\text{sinc}\\left(\\frac{\\pi t}{3}\\right)\\text{sinc}\\left(\\frac{\\pi t}{5}\\right)\\text{sinc}\\left(\\frac{\\pi t}{7}\\right)\\text{sinc}\\left(\\frac{\\pi t}{9}\\right)\\text{sinc}\\left(\\frac{\\pi t}{11}\\right)\\text{sinc}\\left(\\frac{\\pi t}{13}\\right)\\,dt=1"),
            TexMobject("\\displaystyle \\int\\limits_{-\\infty}^{\\infty}\\text{sinc}(\\pi t)\\text{sinc}\\left(\\frac{\\pi t}{3}\\right)\\text{sinc}\\left(\\frac{\\pi t}{5}\\right)\\text{sinc}\\left(\\frac{\\pi t}{7}\\right)\\text{sinc}\\left(\\frac{\\pi t}{9}\\right)\\text{sinc}\\left(\\frac{\\pi t}{11}\\right)\\text{sinc}\\left(\\frac{\\pi t}{13}\\right)\\text{sinc}\\left(\\frac{\\pi t}{15}\\right)\\,dt<1"),
        )
        for txt in borintegs:
            txt.scale_in_place(1 / 3)
        self.play(FadeInFrom(borintegs[0], UP))
        self.wait()
        for i in range(1, len(borintegs)):
            self.play(FadeOutAndShift(borintegs[i - 1], DOWN), FadeInFrom(borintegs[i], UP))
            self.wait()
        self.wait(5)

class SclaedSincFunctions(GraphScene, MovingCameraScene):
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
        self.camera_frame.shift(1 * UP)
        self.play(Write(nplane), ShowCreation(planegrid))
        self.wait()
        def boxcreator(height):
            return VGroup(
                Line(10 * LEFT, LEFT / 2 / height),
                Line(LEFT / 2 / height, LEFT / 2 / height + height * UP),
                Line(LEFT / 2 / height + height * UP, RIGHT / 2 / height + height * UP),
                Line(RIGHT / 2 / height + height * UP, RIGHT / 2 / height),
                Line(RIGHT / 2 / height, 10 * RIGHT)
            ).set_color(BLUE)
        textgroup = VGroup(
            TexMobject("\\displaystyle \\mathcal{F} \\left\\{ \\text{sinc} (\\pi t) \\right\\}=\\text{rect}(\\xi)"),
            TexMobject("\\displaystyle \\mathcal{F}\\left\\{\\frac{\\text{sinc}(\\pi t)}{2}\\right\\}=2\\text{rect}(2\\xi)"),
            TexMobject("\\displaystyle \\mathcal{F}\\left\\{\\frac{\\text{sinc}(\\pi t)}{3}\\right\\}=3\\text{rect}(3\\xi)"),
            TexMobject("\\displaystyle \\mathcal{F}\\left\\{\\frac{\\text{sinc}(\\pi t)}{4}\\right\\}=4\\text{rect}(4\\xi)"),
            TexMobject("\\displaystyle \\mathcal{F}\\left\\{\\frac{\\text{sinc}(\\pi t)}{5}\\right\\}=5\\text{rect}(5\\xi)"),
        )
        for txt in textgroup:
            txt.add_background_rectangle()
        textgroup.shift(1.5 * DOWN)
        rects = VGroup()
        for i in range(1, 1 + len(textgroup)):
            rects.add(boxcreator(i))
        tempc, tempt = textgroup[0].copy(), rects[0].copy()
        self.play(Write(tempc), Write(tempt))
        self.wait()
        for i in range(1, len(textgroup)):
            self.play(
                Transform(tempc, textgroup[i].copy()),
                Transform(tempt, rects[i].copy())
            )
            self.wait()
        self.play(
            FadeOut(tempc),
            FadeOut(tempt),
            self.camera_frame.move_to, ORIGIN,
            self.camera_frame.scale, 1 / 2,
        )
        self.wait(5)

class Borwein(GraphScene, MovingCameraScene):
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
        scaleval = 0.5
        nplane = VGroup(DoubleArrow(8 * LEFT, 8 * RIGHT), DoubleArrow(5.5 * DOWN, 5.5 * UP)).set_color(GREY)
        planegrid = VGroup()
        for i in range(-10, 10 + 1):
            planegrid.add(Line(5 * UP + i * RIGHT, 5 * DOWN + i * RIGHT))
        for i in range(-5, 5 + 1):
            planegrid.add(Line(i * UP + 10 * RIGHT, i * UP + 10 * LEFT))
        planegrid.set_color(GREY)
        planegrid.fade(0.875)
        self.camera_frame.scale(scaleval)
        #self.camera_frame.shift(0.5 * UP)
        self.add(nplane, planegrid)
        self.wait()
        def convolvedgrph(obj1, obj2):
            #make sure the x-coordinates of both the objects are the same
            ob1 = [x[1] for x in obj1]
            ob2 = [x[1] for x in obj2]
            l = len(ob2)
            res = signal.convolve(ob1, ob2) / sum(ob2)
            l = (l - 1) // 2
            res = res[l:]
            res = res[:-l]
            res = [[x[0], y, 0] for x, y in zip(obj2, res)]
            return res
        def get_idx(arr, x):
            return np.searchsorted(arr, x)
        def get_from_x_coordinate(obj, x):
            temp = obj.points
            indx = get_idx(temp[:, 0], x)
            return temp[indx]
        def area_creator(obj1, x1, x2):
            resobj = VMobject(fill_color = YELLOW, fill_opacity = 0.5)
            tpts = obj1.points.copy()
            idx1, idx2 = get_idx(tpts[:, 0], x1), get_idx(tpts[:, 0], x2)
            tpts = tpts[idx1: idx2]
            tpts = np.insert(tpts, 0, [[x1, 0, 0]], axis = 0)
            tpts = np.append(tpts, [[x2, 0, 0], [x1, 0, 0]], axis = 0)
            resobj.set_points_as_corners(tpts)
            return resobj
        def uconv_animator(funcpts, height = 1, rtime = 10, colorpref = GREEN):
            eps = funcpts[1][0] - funcpts[0][0]
            xmi, xma = funcpts[0][0], funcpts[-1][0]
            unitbox = [np.array([x, height if -1 / 2 / height <= x <= 1 / 2 / height else 0, 0]) for x in np.arange(xmi, xma + eps, eps)]
            func_g = VMobject().set_points_as_corners(funcpts)
            movsq = VGroup(
            Line(LEFT / 2 / height, LEFT / 2 / height + height * UP),
            Line(LEFT / 2 / height + height * UP, RIGHT / 2 / height + height * UP),
            Line(RIGHT / 2 / height + height * UP, RIGHT / 2 / height),
            DashedLine(ORIGIN, height * UP),
            )
            convpts = convolvedgrph(funcpts, unitbox)
            convpts_g = VMobject().set_points_as_corners(convpts)
            movsq.shift(xmi * RIGHT)
            func_area = area_creator(func_g.copy().stretch_about_point(height, 1, ORIGIN), movsq.get_critical_point(DOWN)[0] - 1 / 2 / height, movsq.get_critical_point(DOWN)[0] + 1 / 2 / height)
            func_path = VMobject(color = colorpref)
            func_dot = Dot(movsq[-1].get_start(), color = colorpref, radius = DEFAULT_DOT_RADIUS / 3)
            func_path.set_points_as_corners([func_dot.get_center(), func_dot.get_center() + 0.001 * UP])
            func_path_grp = VGroup(func_path, func_dot, convpts_g.fade(1))
            func_grp = VGroup(movsq, func_area, func_g.fade(1))
            def conv_updater(obj):
                tm, ta, tc = obj
                curr = tm.get_critical_point(DOWN)[0]
                ta.become(area_creator(tc.copy().stretch_about_point(height, 1, ORIGIN), curr - 1 / 2 / height, curr + 1 / 2 / height))
            def path_updater(obj):
                pat, pd, pc = obj
                curr = get_from_x_coordinate(pc, movsq[-1].get_start()[0])
                pd.move_to(curr)
                old_path = pat.copy()
                old_path.append_vectorized_mobject(Line(old_path.points[-1], curr))
                pat.become(old_path)
            func_grp.add_updater(conv_updater)
            func_path_grp.add_updater(path_updater)
            self.add(func_grp, func_path_grp)
            self.play(movsq.shift, (xma - xmi) * RIGHT, run_time = rtime, rate_func = linear)
            func_grp.clear_updaters()
            func_path_grp.clear_updaters()
            self.remove(func_grp, func_path_grp)
            return convpts
        fougrp = VGroup(
            TexMobject("\\displaystyle \\mathcal{F}\\left\\{\\text{sinc}(\\pi t)\\right\\}").scale_in_place(scaleval).add_background_rectangle(),
            TexMobject("\\displaystyle \\mathcal{F}\\left\\{\\text{sinc}\\left(\\pi t \\right)\\text{sinc}\\left(\\frac{\\pi t}{3}\\right)\\right\\}").scale_in_place(scaleval).add_background_rectangle(),
            TexMobject("\\displaystyle \\mathcal{F}\\left\\{\\text{sinc}\\left(\\pi t \\right)\\text{sinc}\\left(\\frac{\\pi t}{3}\\right)\\text{sinc}\\left(\\frac{\\pi t}{5}\\right)\\right\\}").scale_in_place(scaleval).add_background_rectangle(),
            TexMobject("\\displaystyle \\mathcal{F}\\left\\{\\text{sinc}\\left(\\pi t \\right)\\text{sinc}\\left(\\frac{\\pi t}{3}\\right)\\text{sinc}\\left(\\frac{\\pi t}{5}\\right)\\text{sinc}\\left(\\frac{\\pi t}{7}\\right)\\right\\}").scale_in_place(scaleval).add_background_rectangle(),
            TexMobject("\\displaystyle \\mathcal{F}\\left\\{\\text{sinc}\\left(\\pi t \\right)\\text{sinc}\\left(\\frac{\\pi t}{3}\\right)\\text{sinc}\\left(\\frac{\\pi t}{5}\\right)\\text{sinc}\\left(\\frac{\\pi t}{7}\\right)\\text{sinc}\\left(\\frac{\\pi t}{9}\\right)\\right\\}").scale_in_place(scaleval).add_background_rectangle(),
            TexMobject("\\displaystyle \\mathcal{F}\\left\\{\\text{sinc}\\left(\\pi t \\right)\\text{sinc}\\left(\\frac{\\pi t}{3}\\right)\\text{sinc}\\left(\\frac{\\pi t}{5}\\right)\\text{sinc}\\left(\\frac{\\pi t}{7}\\right)\\text{sinc}\\left(\\frac{\\pi t}{9}\\right)\\text{sinc}\\left(\\frac{\\pi t}{11}\\right)\\right\\}").scale_in_place(0.4375).add_background_rectangle(),
            TexMobject("\\displaystyle \\mathcal{F}\\left\\{\\text{sinc}\\left(\\pi t \\right)\\text{sinc}\\left(\\frac{\\pi t}{3}\\right)\\text{sinc}\\left(\\frac{\\pi t}{5}\\right)\\text{sinc}\\left(\\frac{\\pi t}{7}\\right)\\text{sinc}\\left(\\frac{\\pi t}{9}\\right)\\text{sinc}\\left(\\frac{\\pi t}{11}\\right)\\text{sinc}\\left(\\frac{\\pi t}{13}\\right)\\right\\}").scale_in_place(0.375).add_background_rectangle(),
            TexMobject("\\displaystyle \\mathcal{F}\\left\\{\\text{sinc}\\left(\\pi t \\right)\\text{sinc}\\left(\\frac{\\pi t}{3}\\right)\\text{sinc}\\left(\\frac{\\pi t}{5}\\right)\\text{sinc}\\left(\\frac{\\pi t}{7}\\right)\\text{sinc}\\left(\\frac{\\pi t}{9}\\right)\\text{sinc}\\left(\\frac{\\pi t}{11}\\right)\\text{sinc}\\left(\\frac{\\pi t}{13}\\right)\\text{sinc}\\left(\\frac{\\pi t}{15}\\right)\\right\\}").scale_in_place(0.3125).add_background_rectangle(),
        )
        fougrp.shift(0.5 * DOWN)
        wid = TextMobject("Width of top edge is 1").scale(scaleval).shift(1.5 * DOWN).add_background_rectangle()
        cgroup = VGroup()
        color_list = [BLUE, GREEN, TEAL, GOLD, RED, PURPLE, ORANGE, PINK, MAROON]
        pts = [np.array([x, 1 if -1 / 2 <= x <= 1 / 2 else 0, 0]) for x in np.arange(-4, 4 + 1 / 75, 1 / 75)]
        grph = VGroup(
            Line(5 * LEFT, LEFT / 2),
            Line(LEFT / 2, LEFT / 2 + UP),
            Line(LEFT / 2 + UP, RIGHT / 2 + UP),
            Line(RIGHT / 2 + UP, RIGHT / 2),
            Line(RIGHT / 2, RIGHT / 2 + 5 * RIGHT),
        )
        grph.set_color(BLUE)
        cgroup.add(grph)
        temptext = fougrp[1].copy()
        widgp = VGroup(
            TexMobject("1-\\frac{1}{3}>0 \\Rightarrow \\mathcal{F} \\{\\cdot\\} \\Big\\vert_{\\xi=0}=1"),
            TexMobject("1-\\frac{1}{3}>0 \\Rightarrow \\mathcal{F} \\{\\cdot\\} \\Big\\vert_{\\xi=0}=1"),
            TexMobject("1-\\frac{1}{3}-\\frac{1}{5}>0 \\Rightarrow \\mathcal{F}\\{\\cdot \\}\\Big\\vert_{\\xi=0}=1"),
            TexMobject("1-\\frac{1}{3}-\\frac{1}{5}-\\frac{1}{7}>0 \\Rightarrow \\mathcal{F}\\{\\cdot \\}\\Big\\vert_{\\xi=0}=1"),
            TexMobject("1-\\frac{1}{3}-\\frac{1}{5}-\\frac{1}{7}-\\frac{1}{9}>0 \\Rightarrow \\mathcal{F}\\{\\cdot \\}\\Big\\vert_{\\xi=0}=1"),
            TexMobject("1-\\frac{1}{3}-\\frac{1}{5}-\\frac{1}{7}-\\frac{1}{9}-\\frac{1}{11}>0 \\Rightarrow \\mathcal{F}\\{\\cdot \\}\\Big\\vert_{\\xi=0}=1"),
            TexMobject("1-\\frac{1}{3}-\\frac{1}{5}-\\frac{1}{7}-\\frac{1}{9}-\\frac{1}{11}-\\frac{1}{13} \\approx 0.045 >0 \\Rightarrow \\mathcal{F}\\{\\cdot \\}\\Big\\vert_{\\xi=0}=1"),
            TexMobject("1-\\frac{1}{3}-\\frac{1}{5}-\\frac{1}{7}-\\frac{1}{9}-\\frac{1}{11}-\\frac{1}{13} -\\frac{1}{15} < 0 \\Rightarrow \\mathcal{F}\\{\\cdot \\}\\Big\\vert_{\\xi=0}<1"),
        )
        for txt in widgp:
            txt.scale_in_place(scaleval)
            txt.add_background_rectangle()
        widgp.shift(1.5 * DOWN)
        wid = widgp[0].copy()
        i = 1
        self.play(
            Write(grph),
            Write(temptext),
        )
        pts = uconv_animator(pts, 2 * i + 1, colorpref = color_list[i + 1])
        cgroup.add(VMobject(color = color_list[i + 1]).set_points_as_corners(pts))
        self.add(cgroup[-1])
        self.wait()
        self.play(
            FadeOut(cgroup[i - 1]),
            Write(wid)
        )
        self.wait()
        for i in range(2, len(fougrp)):
            self.play(
                Transform(temptext, fougrp[i].copy()),
                wid.fade, 0.875,
            )
            pts = uconv_animator(pts, 2 * i + 1, colorpref = color_list[i + 1])
            cgroup.add(VMobject(color = color_list[i + 1]).set_points_as_corners(pts))
            self.add(cgroup[-1])
            self.wait()
            self.play(
                FadeOut(cgroup[i - 1]),
                Transform(wid, widgp[i].copy())
            )
            self.wait()
        self.wait(5)

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
        def convolvedgrph(obj1, obj2):
            #make sure the x-coordinates of both the objects are the same
            ob1 = [x[1] for x in obj1]
            ob2 = [x[1] for x in obj2]
            l = len(ob2)
            res = signal.convolve(ob1, ob2) / sum(ob2)
            l = (l - 1) // 2
            res = res[l:]
            res = res[:-l]
            res = [[x[0], y, 0] for x, y in zip(obj2, res)]
            return res
        eps = 1 / 75
        res = [np.array([x, 1 if -0.5 <= x <= 0.5 else 0, 0]) for x in np.arange(-4, 4 + eps, eps)]
        temp = [np.array([x, 1 if -1 / 2 <= x <= 1 / 2 else 0, 0]) for x in np.arange(-4, 4 + eps, eps)]
        res_g = VMobject().set_points_as_corners(res)
        c_list = [BLUE, GREEN, TEAL, GOLD, RED, PURPLE, ORANGE, PINK, MAROON]
        self.play(
            self.camera_frame.scale, 1 / 2,
            self.camera_frame.shift, 0.5 * UP,
        )
        self.play(ShowCreation(res_g))
        self.wait()
        for i in range(15):
            res = convolvedgrph(res, temp)
            res_g = VMobject(color = c_list[i % 9]).set_points_smoothly(res)
            self.play(ShowCreation(res_g))
            self.wait()
        t = TextMobject("Convolutions and Fourier Transforms").scale(0.75).shift(2 * UP)
        self.play(Write(t))
        self.wait(5)

if __name__ == "__main__":
    module_name = os.path.basename(__file__)
    #command_A = "manim -p  --video_dir ~/Downloads/  "
    clear_cmd = "cls"
    command_A = "python -m manim "
    #command_B = module_name + " " + "SmoothingFunction" + " -pl"
    #command_B = module_name + " " + "ProblemIntroduction" + " -ps -n 16,17"
    #command_B = module_name + " " + "IntroductionScene" + " -p",
    #command_B = module_name + " -p"
    command_B = module_name + " -a"
    #command_B = module_name + " -al"
    os.system(clear_cmd)
    os.system(command_A + command_B)