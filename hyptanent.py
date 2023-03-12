import os
from manim import *


class ExampleCone(ThreeDScene):
    def construct(self):
        axes = ThreeDAxes()
        cone = Cone(base_radius = 1, direction = -Z_AXIS)
        xdot = Dot3D(RIGHT, color = RED)
        ydot = Dot3D(UP, color = GREEN)
        zdot = Dot3D(OUT, color = YELLOW)
        cone.fade = 1 / 2
        self.set_camera_orientation(phi = 60 * DEGREES, theta = 30 * DEGREES, gamma = 0, zoom = 2)
        self.add(axes, xdot, ydot, zdot)
        b = SVGMobject("/Users/muthuveerappanramalingam/Downloads/Manim/HyperbolaTangent/Bubbles_thought")
        self.play(FadeIn(b))
        self.play(Write(cone))
        self.wait(3)

class SVGTest(Scene):
    def construct(self):
        ques2 = Tex("Looks very hard..")
        ques1 = Tex(".. But can be solved \\\\ in closed form")
        bub = SVGMobject("Bubbles_thought")
        bub.set_color(WHITE).set_style(fill_color = BLACK, fill_opacity = 3 / 4, stroke_color = WHITE, stroke_width = 3).scale(2).stretch(2, 0)
        bub.flip()
        def add_and_resize(bubb, contt, content_scale_factor = 0.75, bub_cent_adj_fac = 1 / 8):
            bub, cont = bubb.copy(), contt
            twid = cont.width
            twid += max(2, MED_LARGE_BUFF)
            tht = cont.height
            tht += 2.5 * LARGE_BUFF
            bub.width, bub.height = twid, tht
            bubcent = bub.get_center() + bub_cent_adj_fac * bub.height * UP
            swid = content_scale_factor * bub.width
            if cont.width > swid:
                cont.width = swid
            cont.shift(bubcent - cont.get_center())
            return VGroup(bub, cont)
        ques1grp = add_and_resize(bub.copy(), ques1)
        ques2grp = add_and_resize(bub.copy().flip(), ques2)
        ques1grp.to_corner(DR)
        ques2grp.to_corner(DL)

        self.play(
            Write(ques2grp)
        )
        self.play(Write(ques1grp))
        self.wait()
        self.wait(3)

class Introduction(MovingCameraScene):
    def construct(self):
        self.camera.frame.shift(6 * RIGHT  + 3 * UP )
        eqn = MathTex("f(x)=\\frac{1}{x}").scale(6 / 4).shift(8 * RIGHT + 5 * UP)
        tics = VGroup()
        singletick = Line(DOWN / 8, UP / 8, color = GREY)
        for i in range(2, 1 + 12, 2):
            tics.add(singletick.copy().shift(i * RIGHT))
        for i in range(2, 1 + 6, 2):
            tics.add(singletick.copy().rotate_about_origin(90 * DEGREES).shift(i * UP))
        xyaxis = VGroup(
            Line(DOWN / 2, 6 * UP, color = GREY),
            Line(LEFT / 2, 12 * RIGHT, color = GREY),
            tics
        ).fade(1 / 8)
        self.add(xyaxis)
        self.play(Write(eqn), run_time = 2)
        self.wait()

        hypconst = 8
        hypcurve = FunctionGraph(
            lambda t: hypconst / t,
            x_range = [hypconst / 6, 12],
            color = BLUE
        )
        self.play(Write(hypcurve), run_time = 2)
        self.wait()

        x = ValueTracker(hypconst / 6)
        def create_rect(x0):
            return Rectangle(width = x0, height = hypconst / x0, stroke_width = 1 / 2, color = BLUE, fill_opacity = 3 / 4).shift(x0 * RIGHT / 2 + hypconst / x0 * UP / 2)
        firstrec = create_rect(x.get_value())
        self.play(FadeIn(firstrec))
        self.wait()
        self.add(firstrec)
        firstrec.add_updater(lambda obj: obj.become(create_rect(x.get_value())))

        self.play(
            x.animate.set_value(12),
            run_time = 7,
            rate_func = linear
        )
        self.wait()

        val1, val2 = 4.25, 4.75
        self.play(
            x.animate.set_value(val1),
            run_time = 7 * (12 - val1) / (12 - hypconst / 6),
            rate_func = linear
        )
        self.wait()

        xlabels = VGroup(
            MathTex("x"),
            MathTex("\\frac{1}{x}")
        )
        rectbraces = VGroup(
            Brace(firstrec, UP),
            Brace(firstrec, RIGHT)
        )
        xlabels[0].next_to(rectbraces[0], UP)
        xlabels[1].next_to(rectbraces[1], RIGHT)

        arealabel = MathTex("\\text{Area}=1\\text{ sq. units}")
        arealabel.move_to(firstrec.get_center())
        #arealabel.move_to(6 * RIGHT + 4 * UP)
        #areaarrow = Arrow(firstrec.get_center(), arealabel.get_bottom())
        arealabel.set_background_stroke(width = 3)
        self.play(
            Write(arealabel, run_time = 2),
            #Write(areaarrow)
        )
        self.wait()

        self.play(
            Write(xlabels),
            Write(rectbraces),
            run_time = 3
        )
        self.wait()

        secondrec = create_rect(val2)
        self.bring_to_back(secondrec)
        self.play(
            Write(secondrec),
            xlabels[1].animate.shift((val2 - val1) * RIGHT),
            rectbraces[1].animate.shift((val2 - val1) * RIGHT),
        )
        self.wait()

        diffrects = VGroup(
            Rectangle(width = val1, height = hypconst / val1 - hypconst / val2, stroke_width = 1 / 2, color = RED, fill_opacity = 3 / 4),
            Rectangle(width = val2 - val1, height = hypconst / val2, stroke_width = 1 / 2, color = GREEN, fill_opacity = 3 / 4)
        )
        diffrects[0].move_to(firstrec.get_critical_point(UL) + diffrects[0].get_center() - diffrects[0].get_critical_point(UL))
        diffrects[1].move_to(secondrec.get_critical_point(UR) + diffrects[1].get_center() - diffrects[1].get_critical_point(UR))
        self.play(FadeIn(diffrects))
        self.wait()

        dlines = VGroup(
            DashedLine(diffrects[0].get_critical_point(DR), diffrects[0].get_critical_point(DR) + 1.5 * UP, color = GREEN),
            DashedLine(diffrects[1].get_critical_point(UL), diffrects[0].get_critical_point(DR) + 2 * RIGHT, color = RED),
        )
        dlines.add(dlines[0].copy().shift(RIGHT * (val2 - val1)))
        dlines.add(dlines[1].copy().shift(UP * (hypconst / val1 - hypconst / val2)))
        self.play(Write(dlines), run_time = 2)
        self.wait()

        infbraces = VGroup(
            Brace(Line(dlines[0].get_end(), dlines[2].get_end()), UP, buff = 0),
            Brace(Line(dlines[1].get_end(), dlines[3].get_end()), RIGHT, buff = 0)
        )
        inflabels = VGroup(
            MathTex("dx").set_background_stroke(width = 1, color = WHITE).next_to(infbraces[0], UP),
            #MathTex("d\\left(\\frac{1}{x}\\right)").next_to(infbraces[1], RIGHT)
            MathTex("d\\left(x^{-1}\\right)").set_background_stroke(width = 1, color = WHITE).next_to(infbraces[1], RIGHT)
        )
        inflabels[0].set_color(GREEN)
        inflabels[1].set_color(RED)
        self.play(
            Write(infbraces),
            Write(inflabels),
            run_time = 2
        )
        self.wait()

        comm2 = Tex("Area should always \\\\ be 1 units...")
        comm1 = Tex("... Net change in area \\\\ should be zero")
        bub = SVGMobject("Bubbles_thought")
        bub.set_color(WHITE).set_style(fill_color = BLACK, fill_opacity = 3 / 4, stroke_color = WHITE, stroke_width = 3).scale(2).stretch(2, 0)
        bub.flip()
        def add_and_resize(bubb, contt, content_scale_factor = 0.75, bub_cent_adj_fac = 1 / 8):
            bub, cont = bubb.copy(), contt
            twid = cont.width
            twid += max(2, MED_LARGE_BUFF)
            tht = cont.height
            tht += 2.5 * LARGE_BUFF
            bub.width, bub.height = twid, tht
            bubcent = bub.get_center() + bub_cent_adj_fac * bub.height * UP
            swid = content_scale_factor * bub.width
            if cont.width > swid:
                cont.width = swid
            cont.shift(bubcent - cont.get_center())
            return VGroup(bub, cont)
        ques1grp = add_and_resize(bub.copy(), comm1)
        ques2grp = add_and_resize(bub.copy().flip(), comm2)
        ques1grp.to_corner(DR).shift(6 * RIGHT  + 3 * UP)
        ques2grp.to_corner(DL).shift(6 * RIGHT  + 3 * UP)

        self.play(Write(ques2grp))
        self.play(Write(ques1grp))
        self.wait()

        self.play(
            FadeOut(ques1grp, shift = DOWN),
            FadeOut(ques2grp, shift = DOWN),
            eqn.animate.shift(1 * UP + 4 * LEFT)
        )
        self.wait()

        diffrectscopy = diffrects.copy()
        eqsign = MathTex("=")
        eqsign.move_to(11 * RIGHT + 6 * UP)
        diffrectscopy[0].next_to(eqsign, LEFT)
        diffrectscopy[1].next_to(eqsign, RIGHT)
        diffrectscopy[0].save_state()
        diffrectscopy[1].save_state()
        diffrectscopy[0].move_to(diffrects[0])
        diffrectscopy[1].move_to(diffrects[1])
        self.play(
            Write(eqsign),
            Restore(diffrectscopy[0], path_arc = 135 * DEGREES),
            Restore(diffrectscopy[1], path_arc = -135 * DEGREES),
            run_time = 3
        )
        self.wait()

        eqcopy = eqsign.copy().shift(DL + DOWN)
        inflablescopy = inflabels.copy()
        inflablescopy.add(
            MathTex("-"),
            MathTex("\\cdot", "x"),
            MathTex("\\cdot", "\\frac{1}{x}")
        )
        inflablescopy[3].next_to(eqcopy, LEFT)
        inflablescopy[1].next_to(inflablescopy[3], LEFT, buff = 1 / 8)
        inflablescopy[2].next_to(inflablescopy[1], LEFT, buff = 1 / 8)
        inflablescopy[0].next_to(eqcopy, RIGHT)
        inflablescopy[4].next_to(inflablescopy[0], RIGHT, buff = 1 / 8)
        self.play(Write(inflablescopy, run_time = 3), Write(eqcopy))
        self.wait()

        finaleqn = VGroup(
            MathTex("d(x^{-1})", "\\over", "dx"),
            MathTex("=", "-"),
            MathTex("1", "\\over", "x^2"),
        )
        #finaleqn = MathTex("\\text{Area}", "(", "\\text{region}", ")", "\\over", "\\text{Area}(", "\\text{Circle}", ")")
        finaleqn[0][0].set_color(RED)
        finaleqn[0][-1].set_color(GREEN)
        finaleqn[0][0].set_background_stroke(width = 1, color = WHITE)
        finaleqn[0][-1].set_background_stroke(width = 1, color = WHITE)
        finaleqn.arrange()
        finaleqn.move_to(eqcopy.get_center() + 1.5 * DOWN + finaleqn.get_center() - finaleqn[1][0].get_center())
        self.play(inflabels[-1].animate.shift(1.5 * LEFT + 9 * UP / 16))
        self.wait()

        self.play(TransformMatchingTex(Group(eqcopy, inflablescopy).copy(), finaleqn), run_time = 2)
        firstrec.clear_updaters()
        self.wait()

class HereIsAHyperbola(MovingCameraScene):
    def construct(self):
        self.camera.frame.shift(6 * RIGHT  + 3 * UP)
        tics = VGroup()
        singletick = Line(DOWN / 8, UP / 8, color = GREY)
        for i in range(2, 1 + 12, 2):
            tics.add(singletick.copy().shift(i * RIGHT))
        for i in range(2, 1 + 6, 2):
            tics.add(singletick.copy().rotate_about_origin(90 * DEGREES).shift(i * UP))
        xyaxis = VGroup(
            Line(DOWN / 2, 6 * UP, color = GREY),
            Line(LEFT / 2, 12 * RIGHT, color = GREY),
            tics
        ).fade(1 / 8)
        self.wait()
        self.play(Write(xyaxis))
        self.wait()

        hypconst = 8
        hypcurve = FunctionGraph(
            lambda t: hypconst / t,
            x_range = [hypconst / 6, 12],
            color = BLUE
        )
        self.play(Write(hypcurve), run_time = 2)
        self.wait()

        asymplabel = MathTex("\\text{Asymptotes}").move_to(7 * RIGHT + 4 * UP)
        asymplabeluline = SurroundingRectangle(asymplabel)
        asymparrows = VGroup(
            Arrow(asymplabeluline.get_left(), 3 * UP),
            Arrow(asymplabeluline.get_bottom(), 10 * RIGHT),
        )
        self.play(Write(asymplabel), run_time = 2)
        self.play(Write(asymplabeluline), run_time = 2)
        self.wait()
        self.play(Write(asymparrows), run_time = 2)
        self.wait()
        self.play(LaggedStart(FadeOut(asymplabel), FadeOut(asymplabeluline), FadeOut(asymparrows)))
        self.wait()

        def drawtangent(x0):
            templ = Line(2 * x0 * RIGHT, 2 * UP * hypconst / x0, color = interpolate_color(YELLOW, GREEN, x0 / 5))
            tempgp = VGroup(templ)
            tempgp.add(
                Dot(templ.get_start(), color = PINK, stroke_width = 3, stroke_color = WHITE),
                Dot(templ.get_end(), color = PINK, stroke_width = 3, stroke_color = WHITE),
                Dot(templ.get_center(), color = PINK, stroke_width = 3, stroke_color = WHITE)
            )
            return tempgp
        
        x = ValueTracker(2.5)
        tangentgroup = drawtangent(x.get_value())
        self.play(
            Write(tangentgroup[-1], run_time = 2),
            Write(tangentgroup[:-1])
        )
        self.wait()
        self.add(tangentgroup)

        tangentgroup.add_updater(lambda obj: obj.become(drawtangent(x.get_value())))
        self.play(
            x.animate.set_value(6 - 1 / 16),
            #rate_func = linear,
            run_time = 6
        )
        self.wait()
        apollpre = VGroup(
            Tex("Tangent segment of a hyperbola is ", "bisected", "\\\\ at the tangency point."),
            Tex("- ", "Apollonius", " of Perga")
        )
        apollpre[0][1].set_color(YELLOW)
        apollpre[1][1].set_color(TEAL)
        apollpre.arrange(DOWN)
        apollpre.move_to(8 * RIGHT + 5 * UP)
        apollpre[1].align_to(apollpre[0], RIGHT)
        apollpre.add_background_rectangle()
        self.play(
            Write(apollpre, run_time = 6),
            x.animate.set_value(3.25),
            #rate_func = linear,
            run_time = 4
        )
        self.wait()
        circ = Circle(radius = tangentgroup[0].get_length() / 2).shift(tangentgroup[-1].get_center())
        self.bring_to_back(circ)
        self.play(Write(circ))
        self.wait()
        orig = Dot(color = YELLOW)
        self.play(FadeIn(orig, shift = OUT))
        self.wait()
        labels = VGroup(MathTex("O"), MathTex("P"), MathTex("A"), MathTex("B"))
        labels[0].next_to(orig, direction = UR, buff = 1 / 16)
        labels[1].add_updater(lambda obj: obj.next_to(tangentgroup[-1], direction = UR, buff = 1 / 16))
        labels[2].add_updater(lambda obj: obj.next_to(tangentgroup[1], direction = DOWN, buff = 1 / 8))
        labels[3].add_updater(lambda obj: obj.next_to(tangentgroup[2], direction = LEFT, buff = 1 / 8))
        self.play(LaggedStartMap(Write, labels), run_time = 3)
        self.add(labels[1:])
        self.wait()

        dashradius = DashedLine(ORIGIN, tangentgroup[-1].get_center(), color = GOLD, dash_length = 1 / 8, dashed_ratio = 1 / 2)
        self.play(Write(dashradius))
        self.wait()

        angles = VGroup(
            Angle(Line(ORIGIN, RIGHT), Line(ORIGIN, tangentgroup[-1].get_center()), radius = 1),
            Angle(Line(tangentgroup[1].get_center(), ORIGIN), Line(tangentgroup[1].get_center(), tangentgroup[-1].get_center()), other_angle = True, radius = 1),
        )
        self.play(LaggedStartMap(Write, angles))
        self.wait()
        perpgp = VGroup(
            Dot(3.25 * RIGHT, stroke_width = 3, stroke_color = PINK),
            DashedLine(tangentgroup[-1].get_center(), 3.25 * RIGHT)
        )
        self.play(LaggedStartMap(Write, perpgp), run_time = 2)
        self.wait()

        slopeeqn = VGroup(Tex("tangent's slope ", "="), MathTex("\\text{rise}", "\\over", "\\text{run}"))
        slopeeqn.set_background_stroke(width = 3)
        slopeeqn.arrange(RIGHT)
        slopeeqn.move_to(8 * RIGHT + 2.5 * UP)
        self.play(Write(slopeeqn), run_time = 3)
        self.wait()

        braces = VGroup(
            Brace(Line(ORIGIN, perpgp[0].get_center()), direction = DOWN),
            Brace(Line(perpgp[0].get_center(), tangentgroup[-1].get_center()), RIGHT)
        )
        bracelabels = VGroup(MathTex("x").next_to(braces[0], DOWN), MathTex("\\frac{1}{x}").next_to(braces[1], RIGHT))
        self.play(
            LaggedStartMap(FadeIn, braces, run_time = 2),
            LaggedStartMap(Write, bracelabels, run_time = 2),
        )
        self.wait()

        self.play(
            Rotate(braces[0], PI, about_point = perpgp[0].get_center(), axis = UP),
            bracelabels[0].animate.shift(3.25 * RIGHT),
            run_time = 3
        )
        self.wait()

        fraceqn = VGroup(MathTex("="), MathTex("1/x", "\\over", "-x")).arrange(RIGHT).next_to(slopeeqn, RIGHT)
        #fraceqn = MathTex("=","{1/x}", "\\over", "{-x}").next_to(slopeeqn, RIGHT)
        self.play(TransformMatchingTex(bracelabels.copy(), fraceqn), run_time = 3)
        self.wait()
        self.play(
            FadeOut(circ),
            AnimationGroup(
                FadeOut(braces),
                FadeOut(bracelabels),
                FadeOut(dashradius),
                FadeOut(tangentgroup),
                FadeOut(angles),
                FadeOut(slopeeqn),
                FadeOut(fraceqn),
                FadeOut(perpgp),
                FadeOut(labels),
                FadeOut(orig),
                FadeOut(xyaxis),
                FadeOut(hypcurve),
                lag_ratio = 0.5,
                run_time = 5
            )
        )
        self.wait()

        preplist = VGroup(
            Tex("Preposition 1"),
            Tex("Preposition 2"),
            Tex("Preposition 3"),
            MathTex("\\vdots"),
            Tex("Preposition ..."),
        )
        preplist.arrange(DOWN, buff = 1)
        preplist.align_to(preplist[0].get_critical_point(LEFT), LEFT)
        preplist.set_color(YELLOW)
        preplist[-2].set_color(WHITE)
        preplist.move_to(-0.5 * RIGHT + 6.5 * UP + preplist.get_center() - preplist[0].get_critical_point(UL))
        #self.play(AnimationGroup(*[Write(obj) for obj in preplist], lag_ratio = 0.5), run_time = 8)
        #self.wait()
        contentlist = VGroup(
            Tex("Given an ellipse, circle or a hyperbola ..."),
            Tex("For any ellipse, circle or a hyperbola ..."),
            Tex("Tangent segment of a hyperbola is bisected at the tangency point")
        ).scale(15 / 16)
        contentlist.arrange(DOWN, buff = 1)
        for obj in contentlist[1:]:
            obj.align_to(contentlist[0].get_critical_point(LEFT), LEFT)
        contentlist.move_to(0 * RIGHT + 6 * UP + DOWN / 8 + contentlist.get_center() - contentlist[0].get_critical_point(UL))
        contentlist[-1].shift(3 * DOWN)
        contentlist[2][0][28:36].set_color(BLUE)
        self.play(
            FadeOut(apollpre[2]),
            LaggedStart(
                AnimationGroup(*[Write(obj) for obj in preplist], lag_ratio = 0.5, run_time = 6),
                AnimationGroup(
                    Write(contentlist[:-1]),
                    run_time = 4
                ),
            ),
            Transform(apollpre[1], contentlist[-1], run_time = 2),
        )
        tangentgroup.clear_updaters()
        labels[1].clear_updaters()
        labels[2].clear_updaters()
        labels[3].clear_updaters()
        self.wait()

class ThirdDimension(ThreeDScene):
    def construct(self):
        axes = ThreeDAxes()
        cone = Cone(base_radius = 3, height = 3, direction = -Z_AXIS, fill_color = BLUE, fill_opacity = 7 / 8)
        xdot = Dot3D(RIGHT, color = RED)
        ydot = Dot3D(UP, color = GREEN)
        zdot = Dot3D(OUT, color = YELLOW)
        self.set_camera_orientation(phi = 90 * DEGREES, theta = 0 * DEGREES, gamma = 0, zoom = 3 / 2, frame_center = 0 * OUT)
        self.add(axes, xdot, ydot, zdot)

        hypcurve = FunctionGraph(
            lambda t: np.sqrt(1 + t * t),
            x_range = [-3, 3],
            color = BLUE
        )
        hypcurve.rotate_about_origin(90 * DEGREES, RIGHT)
        hypcurve.rotate_about_origin(90 * DEGREES, OUT)
        hypcurve.scale(2 ** 0.5, about_point = ORIGIN)
        hypcurve.rotate_about_origin(-45 * DEGREES, RIGHT)
        self.play(Write(hypcurve), run_time = 3)
        self.wait()
        self.play(
            hypcurve.animate.scale(1 / 2 ** 0.5, about_point = ORIGIN).rotate_about_origin(45 * DEGREES, RIGHT),
            run_time = 3
        )
        self.wait()
        asymps = VGroup(
            Line(0 * RIGHT + 3 * UP + 3 * OUT, 0 * RIGHT - 3 * UP - 3 * OUT, color = GREEN),
            Line(0 * RIGHT + 3 * UP - 3 * OUT, 0 * RIGHT - 3 * UP + 3 * OUT, color = GREEN)
        )
        self.play(
            Write(asymps),
            hypcurve.animate.shift(RIGHT),
            run_time = 3
        )
        self.wait()

        interline = Line(RIGHT - 3 * UP + 2 * OUT, RIGHT + 3 * UP + 1 * OUT)
        asympinters = VGroup(
            Dot3D(RIGHT - 9/5 * UP + 9/5 * OUT),
            Dot3D(RIGHT + 9/7 * UP + 9/7 * OUT)
        )
        hypinters = VGroup(
            Dot3D(RIGHT - 1.419828 * UP + 1.736638 * OUT),
            Dot3D(RIGHT + 0.905542 * UP + 1.349076 * OUT)
        )
        midpt = Dot3D(RIGHT - 0.2571428 * UP + 1.542857 * OUT)
        asympseg = Line(asympinters[0].get_center(), asympinters[1].get_center(), color = ORANGE)
        hypseg = Line(hypinters[0].get_center(), hypinters[1].get_center(), color = GOLD)
        asympinters.set_color(ORANGE)
        hypinters.set_color(GOLD)
        self.play(Write(interline), run_time = 2)
        self.wait()
        self.play(
            FadeIn(asympinters),
            Write(asympseg, run_time = 2)
        )
        self.wait()
        self.play(
            FadeIn(hypinters),
            Write(hypseg, run_time = 2)
        )
        self.wait()
        self.play(FadeIn(midpt), FadeOut(asympinters), FadeOut(asympseg))
        self.wait()

        asympseg.shift(LEFT)
        asympinters.shift(LEFT)
        hypplane = VMobject()
        hypplane.set_points_as_corners([4 * RIGHT - 3 * UP + 2 * OUT, 4 * RIGHT + 3 * UP + 1 * OUT, -4 * RIGHT + 3 * UP + 1 * OUT, -4 * RIGHT - 3 * UP + 2 * OUT])
        xplane = VMobject(stroke_width = 0, fill_color = BLUE, fill_opacity = 1 / 4)
        xplane.set_points_as_corners([RIGHT + 3 * UP + 3 * OUT, RIGHT - 3 * UP + 3 * OUT, RIGHT - 3 * UP - 3 * OUT, RIGHT + 3 * UP - 3 * OUT, RIGHT + 3 * UP + 3 * OUT])

        self.play(Write(cone), run_time = 3)
        self.wait()

        self.move_camera(phi = 60 * DEGREES, theta = 45 * DEGREES, frame_center = 2 * OUT)
        self.begin_ambient_camera_rotation(rate = 1 / 8)
        self.wait(3)
        self.stop_ambient_camera_rotation()

        self.move_camera(phi = 60 * DEGREES, theta = 45 * DEGREES, frame_center = 2 * OUT)
        self.wait()

        self.play(FadeIn(xplane))
        self.wait(3)

class Conclusion(MovingCameraScene):
    def construct(self):
        self.camera.frame.shift(6 * RIGHT  + 3 * UP )
        tics = VGroup()
        singletick = Line(DOWN / 8, UP / 8, color = GREY)
        for i in range(2, 1 + 12, 2):
            tics.add(singletick.copy().shift(i * RIGHT))
        for i in range(2, 1 + 6, 2):
            tics.add(singletick.copy().rotate_about_origin(90 * DEGREES).shift(i * UP))
        xyaxis = VGroup(
            Line(DOWN / 2, 6 * UP, color = GREY),
            Line(LEFT / 2, 12 * RIGHT, color = GREY),
            tics
        ).fade(1 / 8)
        self.add(xyaxis)
        self.wait()

        eqn = MathTex("f(x)=\\frac{1}{x}").scale(6 / 4).shift(8 * RIGHT + 5 * UP)
        self.play(Write(eqn), run_time = 2)
        self.wait()

        hypconst = 2
        hypcurve = FunctionGraph(
            lambda t: hypconst / t,
            x_range = [hypconst / 6, 12],
            color = BLUE
        )
        self.play(Write(hypcurve), run_time = 2)
        self.wait()

        x, y = ValueTracker(6), ValueTracker(3)
        def create_anyline(x0, y0):
            templ = Line(x0 * RIGHT, y0 * UP, color = ORANGE)
            tempgp = VGroup(templ)
            tempgp.add(Dot(templ.get_start(), stroke_width = 3, color = PINK, stroke_color = WHITE))
            tempgp.add(Dot(templ.get_end(), stroke_width = 3, color = PINK, stroke_color = WHITE))
            return tempgp
        def create_hypchord(x0, y0):
            if x0 * y0 < 4 * hypconst:
                return VGroup()
            temp = (x0 * y0 - 4 * hypconst) ** 0.5
            xyratio = (x0 / y0) ** 0.5
            xymult = (x0 * y0) ** 0.5
            resgp = VGroup()
            fpt = (xymult - temp) * RIGHT / 2 * xyratio + (xymult + temp) * UP / 2 / xyratio
            spt = (xymult + temp) * RIGHT / 2 * xyratio + (xymult - temp) * UP / 2 / xyratio
            resgp.add(
                Line(fpt, spt, color = GOLD),
                Dot(fpt, stroke_width = 3, color = PINK, stroke_color = WHITE),
                Dot(spt, stroke_width = 3, color = PINK, stroke_color = WHITE),
            )
            return resgp
        def create_midpt(x0, y0):
            return Dot(x0 / 2 * RIGHT + y0 / 2 * UP, color = YELLOW, stroke_width = 3, stroke_color = WHITE)
        def create_tangent(x0, y0):
            k = (x0 * y0 / 4 / hypconst) ** 0.5
            x0 /= k
            y0 /= k
            return VGroup(create_anyline(x0, y0), create_midpt(x0, y0))
        aline = create_anyline(x.get_value(), y.get_value())
        cline = create_hypchord(x.get_value(), y.get_value())
        mpt = create_midpt(x.get_value(), y.get_value())
        self.play(
            AnimationGroup(Write(aline), Write(cline), Write(mpt), lag_ratio = 1),
            run_time = 3
        )
        self.wait()

        self.add(aline, cline, mpt)
        aline.add_updater(lambda obj: obj.become(create_anyline(x.get_value(), y.get_value())))
        cline.add_updater(lambda obj: obj.become(create_hypchord(x.get_value(), y.get_value())))
        mpt.add_updater(lambda obj: obj.become(create_midpt(x.get_value(), y.get_value())))

        self.play(
            x.animate.increment_value(4),
            y.animate.increment_value(-1.5),
            run_time = 4,
            #rate_func = linear
        )
        self.play(
            x.animate.increment_value(-4),
            y.animate.increment_value(4),
            #rate_func = linear,
            run_time = 4
        )
        self.wait()
        self.play(
            x.animate.set_value(6),
            y.animate.set_value(4),
            rate_func = linear,
            run_time = 4
        )
        self.wait()
        aline.clear_updaters()
        cline.clear_updaters()
        mpt.clear_updaters()

        xval, yval = x.get_value(), y.get_value()
        ctang = create_tangent(xval, yval)
        r = ValueTracker(0)
        dist = Line(ctang[-1].get_center(), mpt.get_center()).get_length()
        origdist = Line(ctang[-1].get_center(), ORIGIN).get_length()
        def create_movingseg(r0):
            x0 = (origdist + (1 - r0) * dist) * xval / (dist + origdist)
            y0 = (origdist + (1 - r0) * dist) * yval / (dist + origdist)
            resgp = VGroup()
            resgp.add(
                create_anyline(x0, y0),
                create_hypchord(x0, y0),
                create_midpt(x0, y0)
            )
            return resgp
        movingseg = create_movingseg(r.get_value())
        self.add(movingseg)
        movingseg.add_updater(lambda obj: obj.become(create_movingseg(r.get_value())))
        self.play(
            r.animate.set_value(1),
            run_time = 4,
            rate_func = linear
        )
        self.wait()
        movingseg.clear_updaters()
        ctang = create_tangent(x.get_value(), y.get_value())
        self.add(ctang)
        self.remove(movingseg)

        self.play(
            self.camera.frame.animate.scale(3 / 4, about_point = ORIGIN),
            eqn.animate.shift(-(8 * RIGHT + 5 * UP) * (1 / 4)),
            AnimationGroup(
                FadeOut(mpt),
                FadeOut(aline),
                FadeOut(cline),
                lag_ratio = 1 / 3
            )
        )
        self.wait()

        def create_tantriangle(obj):
            x0, y0 = obj.get_start(), obj.get_end()
            resgrp = VGroup(
                Polygon(ORIGIN, x0 / 2, x0 / 2 + y0 / 2, stroke_width = 0 / 4, fill_color = ORANGE, fill_opacity = 1 / 2),
                Polygon(x0 / 2, x0, x0 / 2 + y0 / 2, stroke_width = 0 / 4, fill_color = ORANGE, fill_opacity = 1 / 2),
                Polygon(ORIGIN, y0 / 2, x0 / 2 + y0 / 2, stroke_width = 0 / 4, fill_color = ORANGE, fill_opacity = 1 / 2),
                Polygon(y0 / 2, y0, x0 / 2 + y0 / 2, stroke_width = 0 / 4, fill_color = ORANGE, fill_opacity = 1 / 2)
            )
            return resgrp
        tantriangle = create_tantriangle(ctang[0][0])
        self.play(FadeIn(tantriangle))
        self.wait()
        x0, y0 = ctang[0][0].get_start()[0] / 2, ctang[0][0].get_end()[1] / 2
        constrectangle = Rectangle(height = y0, width = x0).shift(x0 * RIGHT / 2 + y0 * UP / 2)
        coord = MathTex("(x,1/x)").scale(8 / 8)
        coord.next_to(ctang[-1], RIGHT)
        self.play(
            Create(constrectangle),
            Write(coord),
            run_time = 3
        )
        self.wait()

        arealabel = MathTex("\\text{Area } = 2(x)(1/x) =2\\text{ sq. units}").scale(8 / 8).move_to(6 * RIGHT + 2.5 * UP)
        self.play(
            Write(arealabel, run_time = 4),
            Rotate(tantriangle[1], PI, axis = UP, about_point = x0 * RIGHT),
            Rotate(tantriangle[-1], PI, RIGHT, y0 * UP),
            run_time = 4
        )
        self.play(
            Rotate(tantriangle[1], -PI, axis = UP, about_point = x0 * RIGHT),
            Rotate(tantriangle[-1], -PI, RIGHT, y0 * UP),
            run_time = 3
        )
        self.wait()

        self.play(
            AnimationGroup(
                Uncreate(arealabel, run_time = 3),
                FadeOut(coord),
                Uncreate(constrectangle),
                FadeOut(tantriangle),
                run_time = 3,
                lag_ratio = 1 / 3
            )
        )
        self.wait()

        self.add(ctang)
        ctang.add_updater(lambda obj: obj.become(create_tangent(x.get_value(), y.get_value())))
        self.play(
            x.animate.increment_value(6),
            run_time = 3
        )
        self.wait()
        ctang.clear_updaters()

        asympline = Line(DL, 7 * UR, color = GREEN)
        mirrorpt = ctang[-1].copy().rotate_about_origin(PI, UR)

        self.play(
            eqn.animate.shift(RIGHT),
            Write(asympline, run_time = 3),
            FadeTransform(ctang[-1].copy(), mirrorpt, path_arc = 90 * DEGREES, run_time = 3)
        )
        self.wait()

        mirrorline = Line(ORIGIN, mirrorpt.get_center()[0] * RIGHT + mirrorpt.get_center()[1] * UP, color = TEAL)
        dline = DashedLine(ctang[-1].get_center(), mirrorpt.get_center(), dash_length = 1 / 8, dashed_ratio = 1 / 2)
        dline.fade(1 / 2)
        mang = Angle(mirrorline, ctang[0][0], quadrant = [1, -1], elbow = True)
        seccmt = Tex("Radial line", " at a given point and ", "tangent", " \\\\ at a symmetric point are perpendicular")
        seccmt[0].set_color(TEAL)
        seccmt[2].set_color(ORANGE)
        seccmt.scale(7 / 8)
        seccmt.set_background_stroke(width = 3)
        seccmt.move_to(6 * RIGHT + 4 * UP)

        self.play(
            FadeOut(eqn, shift = UP),
            Write(seccmt, run_time = 5),
            Create(mirrorline, run_time = 2),
            Create(dline, run_time = 2),
            FadeIn(mang),
            #asympline.animate.fade(1 / 2),
        )
        self.wait()
        self.play(
            FadeOut(mang),
            FadeOut(dline),
        )
        self.wait()

        ocir = Circle(radius = 2)

        cirhyppt = mirrorline.get_length()
        cirhyppt = mirrorpt.copy().scale(4 / cirhyppt / cirhyppt, about_point = ORIGIN)
        cirhyppt = Dot(cirhyppt.get_center(), stroke_width = 3, color = PURPLE, stroke_color = WHITE)
        thicmt = Tex("These two points are inverses w.r.t \\\\ to the ", "red", " circle").scale(7 / 8)
        thicmt[1].set_color(RED)
        thicmt.set_background_stroke(width = 3)
        thicmt.next_to(seccmt, DOWN, buff = 1)
        twoarrows = VGroup(
            Arrow(cirhyppt.get_center(), thicmt.get_critical_point(LEFT)),
            Arrow(mirrorpt.get_center(), thicmt.get_critical_point(LEFT))
        )
        self.play(
            Write(thicmt, run_time = 4),
            Create(twoarrows, run_time = 2),
            FadeIn(cirhyppt, shift = OUT),
            FadeIn(ocir),
        )
        self.wait()

if __name__ == "__main__":
    module_name = os.path.abspath(__file__)
    #command_A = "manim "+ "-pql" + " " + module_name + " " + "Conclusion" + " -n 0"
    #command_A = "manim "+ "-pql" + " " + module_name
    command_A = "manim "+ "-pqp" + " " + module_name
    #output_location = "C:\ManimCE\media"
    #clear_cmd = "cls"
    #command_A = "manim " + module_name + " " + "RealCase" + " " + "-pql -n 42" + " --media_dir " + output_location
    #command_A = "manim " + module_name + " --media_dir " + output_location + " -pqh"
    #command_A = "manim "+ "-p" + " " + module_name + " " + "SphericalTautochrone" + " -n 0"+ " --renderer=opengl"
    #command_A = "manim "+ "-pqh" + " " + module_name + " " + "RandomPointRandomAngle3" + " -n 13,17"
    #command_A = "manim "+ "-pqh" + " " + module_name + " " + "CircumferenceMethodProbability1" + " -n 0,5" + "  --disable_caching"
    #command_A = "manim "+ "-pqh" + " " + module_name + " --disable_caching"
    #command_A = "manim "+ "-pqh" + " " + module_name
    #command_A = "manim "+ "--write_to_movie" + " " + module_name + " " + "SphericalTautochrone" + " -n 0"+ " --renderer=opengl"
    #os.system(clear_cmd)
    os.system(command_A)