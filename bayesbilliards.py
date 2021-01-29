from manimlib.imports import *
import numpy as np
'''
Create a folder named 'assets' inside the Manim folder.
Create three folders named 'raster_images', 'sounds' and 'svg_images' inside the 'assets' folder.
'''

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
        nplane = VGroup(DoubleArrow(8 * LEFT, 8 * RIGHT), DoubleArrow(8 * DOWN, 8 * UP)).set_color(GREY)
        #self.play(Write(nplane))
        #self.wait()
        planegrid = VGroup()
        for i in range(-10, 10 + 1):
            planegrid.add(Line(8 * UP + i * RIGHT, 8 * DOWN + i * RIGHT))
        for i in range(-5, 5 + 1):
            planegrid.add(Line(i * UP + 10 * RIGHT, i * UP + 10 * LEFT))
        planegrid.set_color(GREY)
        planegrid.fade(0.875)
        self.play(ShowCreation(planegrid), run_time = 3)
        self.wait()
        integ = TexMobject("\displaystyle \int\limits_0^1", "\\binom{n}{k}", "x^k", "(1 - x)^{n - k}", "\,dx")
        integcpy = integ.copy()
        #integl = Underline(integcpy)
        self.play(Write(integ))
        self.wait()
        bith = TexMobject("\\displaystyle \\sum_{j=0}^{n-k}\\binom{n-k}{j}", "(-1)^{j}", "x^j").next_to(integ[2], RIGHT)
        self.play(
            Transform(integ[3], bith.copy()),
            integ[-1].next_to, bith, {"direction": RIGHT})
        self.wait()
        adx = TexMobject("x^{j+k}").shift(UP)
        adx.move_to(integ[3][-1][0].get_center() + adx.get_center() - adx[0][0].get_center())
        #self.play(Write(Dot(integ[3][7].get_center())))
        '''self.play(
            Transform(integ[3][-1], adx.copy()),
            FadeOutAndShift(integ[2], adx.get_center() - integ[2].get_center()),
            integ[-1].shift, RIGHT / 2
        )
        self.wait(0.125)
        self.play(
            integ[3].next_to, integ[1], {"direction": RIGHT},
        )
        self.wait(0.125)
        self.play(
            integ[-1].next_to, integ[3], {"direction": RIGHT},
        )
        self.wait(0.125)'''
        self.play(
            Transform(integ[3][-1], adx.copy()),
            FadeOutAndShift(integ[2], adx.get_center() - integ[2].get_center()),
            ApplyMethod(integ[-1].shift, RIGHT / 2),
        )
        self.wait()
        integmod = VGroup(*[integ[i] for i in [0, 1, 3, 4]])
        integmodc = integmod.copy().shift(1.5 * LEFT)
        integ[-1].add_updater(lambda m: m.next_to(integ[3], RIGHT))
        self.play(
            ApplyMethod(integmod.shift, 1.5 * LEFT),
            ApplyMethod(integ[3].next_to, integmodc[1], {"direction": RIGHT}),
        )
        integ[-1].clear_updaters()
        self.wait()
        self.remove(integ[3][-2])
        self.play(
            integ[0].next_to, integ[3][-5], {"direction": 1 * RIGHT / 8},
            integ[3][-1].shift, 0.5 * RIGHT,
            integ[-1].shift, 0.5 * RIGHT,
        )
        self.wait()
        denom = TexMobject("\\over", "j+k+1").next_to(integ[3][0], direction = RIGHT)
        denom.shift(DOWN / 4)
        vgp = VGroup(integ[0], integ[-1], integ[3][-1])
        self.play(
            Transform(vgp, denom),
            integ[3][-5].next_to, denom, {"direction": UP, "buff": 1 / 8},
        )
        self.wait()
        ulbub, urbub = ThoughtBubble(), ThoughtBubble()
        uldot, urdot = Dot().to_corner(DL), Dot().to_corner(DR)
        ultxt = TextMobject("Can we simplify \\\\ this further?")
        urtxt = TextMobject("Does this have \\\\ an intuitive meaning?")
        ulbub.pin_to(uldot)
        urbub.pin_to(urdot)
        ulbub.add_content(ultxt)
        urbub.add_content(urtxt)
        ulbub.resize_to_content()
        urbub.resize_to_content()
        ulgrp = VGroup(ulbub, ultxt)
        urgrp = VGroup(urbub, urtxt)
        self.play(Write(ulgrp))
        self.wait()
        self.play(
            Write(urgrp),
            FadeOut(ulgrp)
        )
        self.wait()
        self.play(FadeOut(urgrp),)
        self.play(
            FadeOut(integ[1]),
            FadeOut(integ[3][0:6]),
            FadeOut(vgp[0]),
            FadeOut(vgp[2])
        )
        self.wait()
        bayes = ImageMobject("Bayes").scale(2.5).shift(UP)
        bayesname = VGroup(
            TextMobject("Thomas Bayes"),
            TextMobject("(1702 - 1761)")
        )
        bayesname.arrange(DOWN)
        bayesname.next_to(bayes, DOWN)
        bayesname.add_background_rectangle()
        self.play(
            Write(bayesname),
            FadeIn(bayes)
        )
        self.wait()
        tmdot = Dot().to_corner(DR)
        tmbub = ThoughtBubble()
        tmtxt = TexMobject("\\mathbb{P}(B|A)=\\frac{\\mathbb{P}(A|B)\\mathbb{P}(B)}{\mathbb{P}(A)}")
        tmbub.pin_to(tmdot)
        tmbub.add_content(tmtxt)
        tmbub.resize_to_content()
        tmgrp = VGroup(tmbub, tmtxt)
        self.play(Write(tmgrp))
        self.wait()
        integcpy.shift(2.5 * UP)
        self.play(
            FadeOut(tmgrp),
            FadeOutAndShift(bayes, RIGHT),
            FadeOutAndShift(bayesname, RIGHT),
            FadeInFrom(integcpy, UP),
            #FadeInFrom(integl, UP),
        )
        self.wait()
        reasons = VGroup(
            TextMobject("$\\bullet$ Solve this integral using probability"),
            TextMobject("$\\bullet$ Shows the result is independent of $k$")
        )
        reasons.arrange(DOWN, buff = 1)
        reasons.align_on_border(LEFT, buff = 2.5)
        reasons.shift(DOWN / 2)
        self.play(Write(reasons[0]))
        self.wait()
        self.play(Write(reasons[1]))
        self.wait()
        self.wait(5)

class BayesMethod(GraphScene, MovingCameraScene):
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
        #self.play(Write(nplane))
        #self.wait()
        planegrid = VGroup()
        for i in range(-10, 10 + 1):
            planegrid.add(Line(8 * UP + i * RIGHT, 8 * DOWN + i * RIGHT))
        for i in range(-5, 5 + 1):
            planegrid.add(Line(i * UP + 10 * RIGHT, i * UP + 10 * LEFT))
        planegrid.set_color(GREY)
        planegrid.fade(0.875)
        self.play(ShowCreation(planegrid), run_time = 3)
        self.wait()
        num_pts = 9
        spt = 4
        dotgrp = VGroup(*[Dot() for i in range(num_pts)])
        dotgrp.add(Dot(color = GREEN))
        self.play(Write(dotgrp))
        self.wait()
        self.play(ApplyMethod(dotgrp.shift, 3 * UP))
        self.wait()
        gap = 10 / num_pts
        pos = 5
        self.play(
            LaggedStart(*[ApplyMethod(dot.move_to, UP + k * gap * RIGHT + pos * LEFT) for k, dot in zip(range(1 + num_pts), dotgrp)]),
            run_time = 2
        )
        self.wait()
        nline = Line(pos * LEFT, pos * RIGHT).shift(DOWN).set_color(BLUE)
        zer = TexMobject("0").next_to(nline.get_left(), LEFT)
        onl = TexMobject("1").next_to(nline.get_right(), RIGHT)
        self.play(*[Write(obj) for obj in [nline, zer, onl]])
        self.wait()
        random.seed(9)
        randlist = [2 * pos * random.random() - pos for k in range(1 + num_pts)]
        self.wait()
        self.play(LaggedStart(*[ApplyMethod(dot.move_to, p * RIGHT + DOWN) for p, dot in zip(randlist, dotgrp)]), run_time = 2)
        self.wait()
        self.play(*[ApplyMethod(obj.shift, UP) for obj in [nline, zer, onl, dotgrp]])
        self.wait()
        for i in range(6, 10):
            random.seed(7 * i)
            temprandlst = [2 * pos * random.random() - pos for k in range(1 + num_pts)]
            tdgrp = dotgrp.copy()
            for p, dot in zip(temprandlst, tdgrp):
                dot.move_to(p * RIGHT)
            self.play(Transform(dotgrp, tdgrp))
            self.wait(0.25)
        self.wait()
        self.play(Flash(dotgrp[-1]))
        xl = TexMobject("x").scale(0.75).next_to(dotgrp[-1].get_center(), DOWN)
        rxl = TexMobject("1-x").scale(0.75).next_to(dotgrp[-1].get_center(), DOWN)
        self.play(Write(xl))
        self.wait()
        lbrace = Brace(Line(nline.get_start(), dotgrp[-1].get_center()), UP)
        rbrace = Brace(Line(dotgrp[-1].get_center(), nline.get_end()), DOWN)
        rxl.next_to(rbrace, DOWN, buff = 0)
        self.play(
            Write(lbrace),
            Write(rbrace),
            ApplyMethod(xl.move_to, lbrace.get_top() + xl.get_center() - xl.get_bottom() + UP / 8),
            Write(rxl),
        )
        self.wait()
        probinteg = TexMobject("\\int\\limits_0^1", "\\binom{n}{k}", "x^k", "(1-x)^{n-k}", "\\,dx").shift(2 * UP)
        for i in [2, 3, 1]:
            self.play(Write(probinteg[i]))
            self.wait()
        self.wait()
        lbrace.add_updater(lambda m: m.become(Brace(Line(nline.get_start(), dotgrp[-1].get_center()), UP)))
        rbrace.add_updater(lambda m: m.become(Brace(Line(dotgrp[-1].get_center(), nline.get_end()), DOWN)))
        rxl.add_updater(lambda m: m.next_to(rbrace, DOWN, buff = 0))
        xl.add_updater(lambda m: m.move_to(lbrace.get_top() + xl.get_center() - xl.get_bottom() + UP / 8))
        self.play(ApplyMethod(dotgrp[-1].shift, 3 * LEFT))
        self.play(
            ApplyMethod(dotgrp[-1].shift, 4 * RIGHT),
            AnimationGroup(Write(probinteg[0]), Write(probinteg[-1]))
        )
        self.play(ApplyMethod(dotgrp[-1].shift, 3 * LEFT))
        for obj in [lbrace, rbrace, rxl, xl]:
            obj.clear_updaters()
        self.wait()
        self.play(*[FadeOut(obj) for obj in [lbrace, rbrace, rxl, xl]])
        self.wait()
        self.play(
            VGroup(nline, dotgrp, zer, onl).shift, UP,
            probinteg.shift, 0.75 * UP,
        )
        self.wait()
        linelen = 4
        scaleval = 0.75
        templateline = VGroup(Line(ORIGIN, linelen * RIGHT).set_color(BLUE))
        for i in range(1, 1 + 1 + num_pts):
            tdot = Dot().scale(scaleval)
            tdot.shift(linelen * i * RIGHT / (2 + num_pts))
            templateline.add(tdot)
        templateline.shift(linelen * LEFT / 2)
        linelist = VGroup()
        for i in range(1 + num_pts):
            tcpy = templateline.copy()
            tcpy[1 + i].set_color(GREEN)
            if 2 * i < 1 + num_pts:
                tcpy.shift(4 * LEFT)
                tcpy.shift(i * DOWN / 2)
            else:
                tcpy.shift(4 * RIGHT)
                tcpy.shift((i - (1 + num_pts) // 2) * DOWN / 2)
            tcpy.shift(0.5 * DOWN)
            linelist.add(tcpy)
        #animlist = [TransformFromCopy(VGroup(dotgrp, nline), obj) for obj in linelist]
        animlist1 = [TransformFromCopy(dotgrp, obj) for obj in linelist[(1 + num_pts) //2 :]]
        animlist2 = [TransformFromCopy(dotgrp, obj) for obj in linelist[: (1 + num_pts) //2]]
        self.play(
            AnimationGroup(*animlist1, lag_ratio = 0.5, run_time = 3),
            AnimationGroup(*animlist2, lag_ratio = 0.5, run_time = 3),
        )
        self.wait()
        oneover = TexMobject("1", "\\over", "n+1").shift(1 * DOWN)
        ordelab = TextMobject("orderings")
        ordelab.next_to(oneover[-1], DOWN)
        self.play(Write(oneover[-1]), Write(ordelab))
        self.wait()
        self.play(
            FadeOut(ordelab),
            Write(oneover[:-1])
        )
        self.wait()
        tdgrp = VGroup()
        for i in range(1, 1 + 1 + num_pts):
            tdgrp.add(Dot().shift(2 * pos * i * RIGHT / (2 + num_pts)))
        tdgrp[(1 + num_pts) // 2].set_color(GREEN)
        tdgrp.shift(5 * LEFT + 2 * UP)
        eqsign = TexMobject("=").shift(RIGHT)
        self.play(
            Write(eqsign),
            ApplyMethod(oneover.next_to, eqsign, RIGHT),
            ApplyMethod(probinteg.move_to, eqsign.get_critical_point(LEFT) + probinteg.get_center() - probinteg.get_critical_point(RIGHT) + LEFT / 4),
            FadeOut(linelist),
            ApplyMethod(VGroup(nline, dotgrp, zer, onl).shift, UP),
            Transform(dotgrp, tdgrp)
        )
        self.wait(5)

class FirstGeneral(GraphScene, MovingCameraScene):
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
        #self.play(Write(nplane))
        #self.wait()
        planegrid = VGroup()
        for i in range(-10, 10 + 1):
            planegrid.add(Line(8 * UP + i * RIGHT, 8 * DOWN + i * RIGHT))
        for i in range(-5, 5 + 1):
            planegrid.add(Line(i * UP + 10 * RIGHT, i * UP + 10 * LEFT))
        planegrid.set_color(GREY)
        planegrid.fade(0.875)
        self.play(ShowCreation(planegrid), run_time = 3)
        self.wait()
        pos = 5
        dotgrp = VGroup()
        num_pts = 9
        dotgrp.add(Line(ORIGIN, 2 * pos * RIGHT).set_color(BLUE))
        for i in range(1, 1 + 1 + num_pts):
            td = Dot().shift(2 * pos * i * RIGHT / (2 + num_pts))
            dotgrp.add(td)
        '''integlists = VGroup(
            TexMobject("\\int\\limits_0^1", "\\binom{a_1+a_2}{a_1,a_2}", "x_1^{a_1}", "(1-x_1)^{a_2}", "\\,dx_1", "=", "\\frac{1}{a_1+a_2+1}"),
            TexMobject("\\int\\limits_0^1\\int\\limits_0^{x_2}", "\\binom{a_1+a_2+a_3}{a_1,a_2,a_3}", "x_1^{a_1}", "(x_2-x_1)^{a_2}", "(1-x_2)^{a_3}", "\\,dx_1dx_2", "=", "\\prod_{k=1}^2\\frac{1}{a_1+a_2+a_3+k}"),
        )'''
        integlists = VGroup(
            TexMobject("\\int\\limits_0^1", "\\binom{n}{a_1,a_2}", "x_1^{a_1}", "(1-x_1)^{a_2}", "\\,dx_1", "=", "\\prod_{k=1}^1\\frac{1}{n+k}"),
            TexMobject("\\int\\limits_0^1\\int\\limits_0^{x_2}", "\\binom{n}{a_1,a_2,a_3}", "x_1^{a_1}", "(x_2-x_1)^{a_2}", "(1-x_2)^{a_3}", "\\,dx_1dx_2", "=", "\\prod_{k=1}^2\\frac{1}{n+k}"),
            TexMobject("\\int\\limits_0^1\\int\\limits_0^{x_3}\\int\\limits_0^{x_2}", "\\binom{n}{a_1,\\cdots,a_4}", "x_1^{a_1}", "(x_2-x_1)^{a_2}", "(x_3-x_2)^{a_3}", "(1-x_3)^{a_4}", "\\,dx_1dx_2dx_3", "=", "\\prod_{k=1}^3\\frac{1}{n+k}"),
            TexMobject("\\int\\limits_0^1\\int\\limits_0^{x_4}\\int\\limits_0^{x_3}\\int\\limits_0^{x_2}", "\\binom{n}{a_1,\\cdots,a_5}", "x_1^{a_1}", "(x_2-x_1)^{a_2}", "(x_3-x_2)^{a_3}", "(x_4-x_3)^{a_4}", "(1-x_4)^{a_5}", "\\,dx_1dx_2dx_3dx_4", "=", "\\prod_{k=1}^4\\frac{1}{n+k}"),
        )
        integlists.arrange(DOWN)
        for txt in integlists:
            txt.scale(0.625)
        for txt in integlists[:2]:
            self.play(Write(txt))
            self.wait()
        dotgrp.move_to(integlists[2].get_center())
        fpos = (1 + 1 + num_pts) // 3
        spos = 2 * (1 + 1 + num_pts) // 3
        dotgrp[fpos].set_color(GREEN)
        dotgrp[spos].set_color(GREEN)
        self.play(GrowFromPoint(dotgrp, UP))
        xlabs = VGroup(
            TexMobject("x_1").scale(0.625).next_to(dotgrp[fpos], DOWN, buff = 1 / 8),
            TexMobject("x_2").scale(0.625).next_to(dotgrp[spos], DOWN, buff = 1 / 8),
        )
        self.play(Write(xlabs))
        self.wait()
        self.play(
            FadeOutAndShiftDown(VGroup(dotgrp, xlabs))
        )
        for txt in integlists[2:]:
            self.play(Write(txt))
            self.wait()
        self.wait(5)

class SecondGeneral(GraphScene, MovingCameraScene):
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
        #self.play(Write(nplane))
        #self.wait()
        planegrid = VGroup()
        for i in range(-10, 10 + 1):
            planegrid.add(Line(8 * UP + i * RIGHT, 8 * DOWN + i * RIGHT))
        for i in range(-5, 5 + 1):
            planegrid.add(Line(i * UP + 10 * RIGHT, i * UP + 10 * LEFT))
        planegrid.set_color(GREY)
        planegrid.fade(0.875)
        #self.play(ShowCreation(planegrid), run_time = 3)
        self.wait()
        pos = 5
        numpts = 9
        random.seed(7)
        xrands = [pos * random.random() for k in range(1 + numpts)]
        yrands = [pos * random.random() for k in range(1 + numpts)]
        dotgrp = VGroup(RegularPolygon(n = 4, start_angle = -PI / 4 - PI / 2).scale(pos / np.sqrt(2)).shift(pos * UR / 2))
        for i in range(1 + numpts):
            dotgrp.add(Dot().shift(xrands[i] * RIGHT + yrands[i] * UP))
        dotgrp.shift(pos * DL / 2)
        dotgrp[-1].set_color(GREEN)
        dotgrp.rotate_about_origin(90 * DEGREES, axis = RIGHT)
        self.play(Write(dotgrp))
        self.play(
            ApplyMethod(dotgrp.rotate_about_origin, -90 * DEGREES, {"axis": RIGHT})
        )
        self.wait()
        lc = dotgrp[0].get_start()
        rc = dotgrp[-1].get_center()
        rect = Rectangle(width = rc[0] - lc[0], height = rc[1] - lc[1], fill_color = GREEN, fill_opacity = 0.25, stroke_opacity = 0.0)
        rect.move_to(dotgrp[0].get_start() + rect.get_center() - rect.get_corner(DL))
        self.play(Write(rect))
        self.wait()
        dartgrp = VGroup(dotgrp, rect)
        self.play(dartgrp.shift, UP)
        probinteg = TexMobject("\\int\\limits_0^1\\int\\limits_0^1", "\\binom{n}{k}", "(xy)^k", "(1-xy)^{n-k}", "\\,dydx").scale(0.875).shift(2.5 * DOWN)
        coords = TexMobject("(x,y)").set_color(GREEN).scale(0.75).next_to(dotgrp[-1].get_center(), RIGHT)
        self.play(Write(coords))
        self.wait()
        for i in [2, 3, 1]:
            self.play(Write(probinteg[i]))
            self.wait()
        self.play(Write(probinteg[0]), Write(probinteg[-1]))
        self.wait()
        self.play(
            ApplyMethod(VGroup(dotgrp, probinteg, coords, rect).shift, 3.5 * LEFT)
        )
        def relordering(l0):
            array = np.array(l0.copy())
            temp = array.argsort()
            ranks = np.empty_like(temp)
            ranks[temp] = np.arange(len(array))
            return ranks
        def twodordering(l0, pos0):
            tdots = l0.copy()
            tpos = pos0
            xord, yord = [], []
            for obj in tdots[1:]:
                tempc = obj.get_center()
                xord += [tempc[0]]
                yord += [tempc[1]]
            xord = relordering(xord)
            yord = relordering(yord)
            tstart = tdots[0].get_start()
            tlen = len(tdots)
            for i in range(1, tlen):
                obj = tdots[i]
                obj.move_to(tstart)
                obj.shift((1 + xord[i - 1]) * pos * RIGHT / tlen + (1 + yord[i - 1]) * pos * UP / tlen)
            return tdots
        def gridcreator(pos, npts):
            res = VGroup()
            for i in range(npts):
                res.add(Line(ORIGIN, pos * RIGHT, stroke_width = 0.5).shift((1 + i) * pos * UP / (1 + npts)))
                res.add(Line(ORIGIN, pos * UP, stroke_width = 0.5).shift((1 + i) * pos * RIGHT / (1 + npts)))
            return res
        mdotgrp = twodordering(dotgrp, pos).shift(7 * RIGHT)
        mgrid = gridcreator(pos, 1 + numpts).move_to(mdotgrp[0].get_center())
        self.play(Write(mdotgrp[0]))
        self.wait()
        self.play(Write(mgrid))
        self.wait()
        anim = [TransformFromCopy(obj1, obj2) for obj1, obj2 in zip(dotgrp[1:], mdotgrp[1:])]
        self.play(FadeOut(coords))
        self.play(
            AnimationGroup(*anim, lag_ratio = 0.05),
            run_time = 3
        )
        self.wait()
        parrow = ArrowTip(color = GREEN).flip().next_to(mdotgrp[-1].get_center(), LEFT, buff = 1 / 16).shift(7 * pos * LEFT / (2 + numpts))
        mlc = mdotgrp[0].get_start()
        mrc = mdotgrp[-1].get_center()
        mrect = Rectangle(width = mrc[0] - mlc[0], height = mrc[1] - mlc[1], fill_color = GREEN, fill_opacity = 0.25, stroke_opacity = 0.0)
        mrect.move_to(mdotgrp[0].get_start() + mrect.get_center() - mrect.get_corner(DL))
        #combeqn = TexMobject("\\frac{1}{n+1}", "\\cdot", "\\sum_{j = k}^n \\frac{1}{j+1}").scale(0.875).shift(2.5 * DOWN + 3.5 * RIGHT)
        combeqn = TexMobject("\\frac{1}{n+1}", "\\cdot", "\\sum_{j = k}^{n} \\frac{1}{j + 1}").scale(0.875).shift(2.5 * DOWN + 3.5 * RIGHT)
        self.play(Write(mrect), Write(parrow))
        self.wait()
        self.play(Write(combeqn[0]))
        self.wait()
        self.play(Write(combeqn[1:]))
        self.wait()
        self.play(FadeOut(parrow), FadeOut(mrect))
        #self.play(Transform(combeqn[-1], combsum))
        eqsign = TexMobject("=").scale(0.875).move_to(2.5 * DOWN + 0.5 * RIGHT)
        self.play(
            Write(eqsign),
            combeqn.next_to, eqsign, RIGHT,
            probinteg.next_to, eqsign, LEFT,
        )
        self.wait()
        self.play(
            FadeOut(VGroup(dotgrp, mdotgrp, eqsign, combeqn, probinteg, rect, mgrid))
        )
        finalstat = TextMobject("Extending this idea to higher dimensional cubes \\\\ with more special darts is straightforward").shift(UP)
        finaleqn = TexMobject("\\int\\limits_0^1\\int\\limits_0^1\\int\\limits_0^1", "\\binom{n}{k}", "(xyz)^k", "(1-xyz)^{n-k}", "\\,dzdydx", "=", "?").shift(DOWN)
        self.play(Write(finalstat))
        self.play(Write(finaleqn))
        self.wait(5)

class AnimationTransform(Scene):
    def construct(self):
        square = Square()
        circle = Circle()
        anno = TextMobject("Transform")
        anno.shift(2 * DOWN)
        self.add(anno)
        self.add(square)
        self.play(Transform(square, circle))
        #square.generate_target()
        #square.target.move_to(2 * UP)
        #self.play(MoveToTarget(square))
        self.play(square.shift, 2 * UP)
        self.wait()
        self.play(circle.shift, UP)
        self.wait(5)

if __name__ == "__main__":
    module_name = os.path.basename(__file__)
    #command_A = "manim -p  --video_dir ~/Downloads/  "
    clear_cmd = "cls"
    command_A = "python -m manim "
    #command_B = module_name + " " + "SecondGeneral" + " -pl -n 21,24"
    #command_B = module_name + " " + "ProblemIntroduction" + " -ps -n 16,17"
    #command_B = module_name + " " + "AnimationTransform" + " -pl"
    command_B = module_name + " -p"
    #command_B = module_name + " -pl"
    #command_B = module_name + " -a"
    #command_B = module_name + " -al"
    os.system(clear_cmd)
    os.system(command_A + command_B)