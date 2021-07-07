from manim import *
import os

from manim.utils.file_ops import add_version_before_extension
from numpy import dot


class InterviewQuestion(Scene):
    def construct(self):
        alice = SVGMobject("/Users/muthuveerappanramalingam/Downloads/ManimInstall/manim_ce/GeometricProbability/assets/alice")
        bob = SVGMobject("/Users/muthuveerappanramalingam/Downloads/ManimInstall/manim_ce/GeometricProbability/assets/bob")
        title = Tex(r"Interview Question")
        utitle = Underline(title)
        titlegrp = VGroup(title, utitle)
        self.play(Write(titlegrp))
        self.play(titlegrp.animate.to_edge(UP))
        self.wait()

        self.play(Write(alice))
        self.play(alice.animate.shift(4.0 * LEFT).scale(1.5))
        self.play(Write(bob))
        self.play(bob.animate.shift(4.0 * RIGHT).scale(1.5))
        self.wait()

        clockrad = 1
        clockgrp = VGroup()
        rimgrp = VGroup()
        cols = [BLUE, GREEN, YELLOW, ORANGE, RED]
        irad, inc = clockrad, 1 / 32
        for col in cols:
            rimgrp.add(Circle(radius = irad, color = col))
            irad += inc
        clockgrp.add(rimgrp)
        hourtick = Line(0.85 * RIGHT, 0.95 * RIGHT)
        tickgrp = VGroup()
        for hr in range(12):
            tickgrp.add(hourtick.copy().rotate_about_origin(hr * 30 * DEGREES))
        clockgrp.add(tickgrp)
        hhand = Line(0.1 * LEFT, 0.4 * RIGHT).set_color(BLUE)
        mhand = Line(0.1 * LEFT, 0.8 * RIGHT).rotate_about_origin(90 * DEGREES).set_color(GREEN)
        clockgrp.add(hhand, mhand)
        clockgrp.scale(1.5)
        hhand.save_state()
        mhand.save_state()
        self.play(Write(clockgrp))

        t = ValueTracker(0)
        mhand.add_updater(lambda m: m.restore().rotate_about_origin(-t.get_value() * 360 * DEGREES))
        hhand.add_updater(lambda m: m.restore().rotate_about_origin(-t.get_value() * 30 * DEGREES))
        self.add(mhand, hhand)
        self.play(
            t.animate.set_value(1),
            run_time = 3,
            rate_func = linear
        )
        mhand.clear_updaters()
        hhand.clear_updaters()
        self.wait()

        ques = MathTex("\\mathbb{P}(\\text{Alice and Bob meet in the library})")
        ques.shift(3 * DOWN)
        self.play(Write(ques))
        self.wait()

        tempaxes = VGroup(Line(0.5 * LEFT, 7 * RIGHT).set_color(GREEN), Line(0.5 * DOWN, 7 * UP).set_color(RED)).shift(3.5 * DL)
        tempaxes.add(Line(3.5 * DOWN + 1.5 * RIGHT, 3.5 * DOWN + 1.5 * RIGHT + DOWN / 8))
        tempaxes.add(Line(3.5 * LEFT + 1.5 * UP, 3.5 * LEFT + 1.5 * UP + LEFT / 8))
        self.play(
            Write(tempaxes),
            FadeOut(clockgrp),
            FadeOut(titlegrp),
            ques.animate.to_corner(UR),
            alice.animate.move_to(3.5 * LEFT + 1.5 * UP).shift(1 * UL + UP / 4).scale(0.75),
            bob.animate.move_to(3.5 * DOWN + 1.5 * RIGHT).shift(1 * UR + UP / 4).scale(0.75),
        )
        self.wait()

        fdot = Dot(5 * 0.75 * RIGHT + 5 * 0.5 * UP).shift(3.5 * DL)
        ftext = MathTex("(\\text{3:45,3:30})").next_to(fdot, UP, buff = 1 / 8)
        ftext[0][1:5].set_color(GREEN)
        ftext[0][6:10].set_color(RED)
        self.play(Write(fdot))
        self.play(Write(ftext))
        self.wait()
        self.play(FadeOut(ftext))
        self.wait()

        anim_list = []
        meetgrp, nmeetgrp = VGroup(), VGroup()
        for k in range(600):
            xord, yord = 5 * np.random.random(), 5 * np.random.random()
            tdot = Dot(xord * RIGHT + yord * UP).shift(3.5 * DL)
            if -5 <= 6 * (xord - yord) <= 5:
                meetgrp.add(tdot)
            else:
                nmeetgrp.add(tdot)
            anim_list += [GrowFromCenter(tdot)]
        self.play(LaggedStart(*anim_list), run_time = 3)
        self.remove(fdot)
        self.wait()
        self.play(
            meetgrp.animate.set_color(BLUE),
            nmeetgrp.animate.set_color(ORANGE)
        )
        self.wait()

        nmeetpoly = VMobject(fill_color = ORANGE, fill_opacity = 1 / 2, stroke_opacity = 0)
        nmeetpoly.set_points_as_corners([5 * RIGHT / 6, 5 * RIGHT, 5 * RIGHT + 25 * UP / 6, 5 * RIGHT / 6])
        nmeetpolycpy = nmeetpoly.copy().rotate_about_origin(180 * DEGREES, axis = UR)
        meetpoly = VMobject(fill_color = BLUE, fill_opacity = 1 / 2, stroke_opacity = 0)
        meetpoly.set_points_as_corners([ORIGIN, 5 * RIGHT / 6, 5 * RIGHT + 25 * UP / 6, 5 * UR, 5 * UR + 5 * LEFT / 6, 5 * UP / 6, ORIGIN])
        meetpolygrp = VGroup(nmeetpoly, nmeetpolycpy, meetpoly).shift(3.5 * DL)

        self.play(Write(meetpolygrp))
        self.wait()
        self.play(
            meetpolygrp.animate.set_style(fill_opacity = 0.75),
            meetgrp.animate.set_opacity(0.75 * 0.5),
            nmeetgrp.animate.set_opacity(0.75 * 0.5),
        )
        self.wait()

        prob = VGroup(
            MathTex("=\\frac{11/36}{11/36+25/36}"),
            MathTex("=\\frac{11}{36}")
        )
        prob[0][0][1:6].set_color(BLUE)
        prob[0][0][7:12].set_color(BLUE)
        prob[0][0][13:18].set_color(ORANGE)
        prob[1][0][1:].set_color(BLUE)
        prob.next_to(ques, DOWN)
        prob[0].to_edge(RIGHT)
        self.play(Write(prob[0]))
        self.wait()

        pholder = ques.copy().move_to(ORIGIN).shift(2.0 * LEFT)
        prob[1].next_to(pholder, RIGHT)
        self.play(
            FadeOut(tempaxes),
            FadeOut(meetpolygrp),
            FadeOut(meetgrp), FadeOut(nmeetgrp),
            FadeOut(alice, shift = LEFT),
            FadeOut(bob, shift = RIGHT),
            ques.animate.move_to(ORIGIN).shift(2.0 * LEFT),
            prob[0].animate.next_to(pholder, RIGHT)
        )
        self.wait()
        self.play(
            FadeOut(prob[0], shift = DOWN),
            FadeIn(prob[1], shift = DOWN)
        )
        self.wait()

        bub = SVGMobject("/Users/muthuveerappanramalingam/Downloads/ManimInstall/manim_ce/GeometricProbability/assets/Bubbles_thought")
        bub.set_color(WHITE).set_style(fill_color = BLACK, fill_opacity = 3 / 4, stroke_color = WHITE, stroke_width = 3).scale_in_place(2).stretch_in_place(2, 0)
        bub.flip()
        nothgeom = Tex("Nothing geometric about \\\\ the problem!")
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
        nothgrp = add_and_resize(bub.copy(), nothgeom)
        nothgrp.to_corner(DR)
        bubf = bub.flip()
        solgeom = Tex("but the solution is \\\\ purely geometric...")
        solgeomgrp = add_and_resize(bubf, solgeom)
        solgeomgrp.to_corner(DL)
        self.play(
            Write(nothgrp),
            ques.animate.shift(2 * UP),
            prob[1].animate.shift(2 * UP)
        )
        self.wait()
        self.play(Write(solgeomgrp))
        self.wait(5)


class GeomProb(Scene):
    def construct(self):
        titletxt = Tex(r"Geometric Probability")
        titleline = Underline(titletxt)
        titlegrp = VGroup(titletxt, titleline)
        self.play(Write(titletxt))
        self.play(Write(titleline))
        self.wait()
        self.play(titlegrp.animate.to_edge(UP))
        self.wait()

        '''result = VGroup(Integer(10, color = GREEN), Tex(r"outter of"), Integer(10)).arrange(RIGHT)
        self.add(result)
        self.play(
            ChangeDecimalToValue(result[0], 0),
            rate_func=lambda t: smooth(t, 3),)'''

        sampgrp = VGroup()
        ccol, scol = GREEN, RED
        rad = 2.5
        for i in range(1000):
            x, y = rad * (2 * np.random.random() - 1), rad * (2 * np.random.random() - 1)
            temp = Dot(x * RIGHT + y * UP)
            if x * x + y * y <= rad * rad:
                temp.set_color(ccol)
            else:
                temp.set_color(scol)
            sampgrp.add(temp)
        cir, sq = Circle(radius = rad, color = ccol), Square(side_length = 2 * rad, color = scol)
        probtxt = MathTex("\\mathbb{P}(\\text{random dot inside a square lies within the circle})").shift(3 * DOWN)
        self.play(Write(cir), Write(sq))
        self.wait()
        self.play(
            AnimationGroup(
                LaggedStart(*[FadeIn(dot, scale = 3) for dot in sampgrp], run_time = 3),
                Write(probtxt),
                lag_ratio = 0.5
            )
        )
        self.wait()
        self.play(VGroup(sampgrp, cir, sq, probtxt).animate.shift(10 * DOWN))
        self.remove(VGroup(sampgrp, cir, sq, probtxt))
        self.wait()

        cirr = Circle(radius = rad, color = BLUE)
        pdot = Dot(rad * DOWN)

        '''t = ValueTracker(0.875)
        def makeqgrp(t):
            qdot = Dot(rad * np.cos(2 * PI * t) * RIGHT + rad * np.sin(2 * PI * t) * UP)
            #qlab.add_background_rectangle()
            qline = Line(pdot.get_center(), qdot.get_center())
            pqnorm = qline.get_length()
            pq = (qdot.get_center() - pdot.get_center()) / (1 if pqnorm <= 0.0000001 else pqnorm)
            qlab = MathTex("Q").scale(0.875).next_to(qdot, pq, buff = 1 / 16)
            #return VGroup(qline, qdot, qlab)
            return VGroup(qdot, qlab, qline)
        def qgrpupd(obj):
            thet = t.get_value()
            obj.become(makeqgrp(thet))
        qgrp = makeqgrp(t.get_value())
        self.play(Write(cirr), Write(pdot), Write(plab), Write(qgrp))
        self.wait()
        self.add(qgrp)
        self.play(t.animate.increment_value(0.25), UpdateFromFunc(qgrp, qgrpupd), run_time = 2)
        self.play(t.animate.increment_value(0.41), UpdateFromFunc(qgrp, qgrpupd), run_time = 2)
        self.play(t.animate.increment_value(0.15), UpdateFromFunc(qgrp, qgrpupd), run_time = 2)'''

        dotgrp = VGroup()
        expectxt = MathTex("\\mathbb{E}(\\text{length of a random chord on a circle})").shift(3 * DOWN)
        for k in range(314):
            t = np.random.random()
            x, y = rad * np.cos(2 * PI * t), rad * np.sin(2 * PI * t)
            #lcol = interpolate_color(GREEN, RED, x / 2 / rad)
            lcol = random_color()
            #x, y = 2 * x - 1, 2 * y - 1
            ddot = Dot(x * RIGHT + y * UP, color = lcol)
            dline = Line(pdot.get_center(), ddot.get_center(), color = lcol)
            dotgrp.add(VGroup(ddot, dline))
        self.play(Write(cirr), Write(pdot))
        self.wait()
        #self.play(ShowIncreasingSubsets(dotgrp))
        self.play(
            AnimationGroup(
                LaggedStart(*[Create(grp) for grp in dotgrp], run_time = 3),
                Write(expectxt),
                lag_ratio = 0.5
            )
        )
        self.wait()
        self.play(FadeOut(cirr), FadeOut(pdot), FadeOut(dotgrp, shift = DOWN), FadeOut(expectxt, shift = DOWN))
        self.wait()

        putn = Tex(r"B1 Putnam 1989").move_to(titletxt)
        putnu = Underline(putn)
        putgrp = VGroup(putn, putnu)
        sq = Square(side_length = 2 * rad, color = BLUE)
        sampgrp = VGroup()
        for k in range(1000):
            x, y = rad * (2 * np.random.random() - 1), rad * (2 * np.random.random() - 1)
            disted = min(rad - x, rad + x, rad - y, rad + y)
            #disted = rad - x
            distcent = np.sqrt(x * x + y * y)
            ddot = Dot(x * RIGHT + y * UP, color = BLUE, radius = 1 / 16)
            if distcent <= disted:
                ddot.set_color(GREEN)
            sampgrp.add(ddot)
        self.play(Transform(titlegrp, putgrp), Write(sq))
        parab = ParametricFunction(lambda t: np.array([t, (t * t - rad * rad) / 2 / rad, 0]), t_range = np.array([-rad, rad]), color = RED)
        parab2 = parab.copy().rotate_about_origin(90 * DEGREES)
        self.wait()
        self.play(LaggedStart(*[FadeIn(dot, scale = 5) for dot in sampgrp], run_time = 3))
        #self.add(parab, parab2)
        probab = VGroup(MathTex("\\mathbb{P}(\\text{random point closer to center than the edges})="), MathTex("\\frac{4\\sqrt{2}-5}{3}"))
        for txt in probab:
            txt.scale(0.875)
        probab.arrange(RIGHT)
        probab.shift(3 * DOWN + LEFT)
        probabn = MathTex("\\frac{1}{12}\\left(4-\\sec^4\\left(\\frac{\\pi}{2n}\\right)\\right)").scale(0.875).next_to(probab[0], RIGHT)
        self.play(Write(probab))
        self.wait()
        self.play(Transform(probab[1], probabn))
        self.wait()
        self.play(
            FadeOut(sq, sampgrp),
            probab.animate.shift(3 * UP),
            FadeOut(putgrp, shift = UP)
        )
        self.wait(3)


class TrickyAndDifficult(Scene):
    def construct(self):
        titletxt = Tex(r"Geometric Probability")
        titleline = Underline(titletxt)
        titlegrp = VGroup(titletxt, titleline)
        titlegrp.to_edge(UP)
        self.add(titlegrp)
        cmar = VMobject(fill_opacity = 1, stroke_color = GREEN, fill_color = GREEN)
        cmar.set_points_as_corners([ORIGIN, 2 * UR + 0.5 * RIGHT, 0.75 * UP, UL, ORIGIN])
        cmar.scale(1 / 8)
        cmargrp = VGroup()
        for k in range(4):
            cmargrp.add(cmar.copy())
        cmargrp.arrange(DOWN, buff = 1)
        txtgrp = VGroup(
            Tex("Simpler to state"),
            Tex("Tricky"),
            Tex("Interesting"),
            Tex("Difficult"),
        )
        cmartxt = VGroup()
        for cm, tx in zip(cmargrp, txtgrp):
            tx.next_to(cm, RIGHT, buff = 1 / 2)
            cmartxt.add(VGroup(cm, tx))
        cmartxt.move_to(ORIGIN)
        self.play(LaggedStart(*[Write(obj) for obj in cmartxt]))
        self.wait()
        self.play(cmartxt.animate.to_edge(LEFT))
        self.wait()

        rad = 2
        mv = 2 * RIGHT
        cir = Circle(radius = rad, color = BLUE).shift(mv)
        def threepts():
            t1, t2, t3 = np.random.random(), np.random.random(), np.random.random()
            t1, t2, t3 = min(t1, t2, t3), t1 + t2 + t3, max(t1, t2, t3)
            t2 -= (t1 + t3)
            arcgrp = VGroup()
            t1, t2, t3 = t1 * 2 * PI, t2 * 2 * PI, t3 * 2 * PI
            arcgrp.add(Arc(radius = rad, start_angle = t1, angle = t2 - t1, color = TEAL))
            arcgrp.add(Arc(radius = rad, start_angle = t2, angle = t3 - t2, color = ORANGE))
            arcgrp.add(Arc(radius = rad, start_angle = t3, angle = 2 * PI + t1 - t3, color = PINK))
            for t in [t1, t2, t3]:
                arcgrp.add(Dot(rad * RIGHT).rotate_about_origin(t))
            return arcgrp.rotate_about_origin(np.random.random() * 2 * PI).shift(mv)
        self.play(Write(cir))
        self.wait()
        tarcs = VGroup(*[threepts() for k in range(15)])
        self.play(Write(tarcs[0]))
        self.wait()
        impdot = Dot(rad * RIGHT, color = BLUE).scale(5 / 4).shift(mv)
        self.bring_to_front(impdot)
        self.play(Write(impdot))
        self.wait()
        for obj in tarcs[1:]:
            self.play(Transform(tarcs[0], obj))
        self.remove(cir)
        self.wait()
        expec = Tex("$\\mathbb{E}$(Length of the arc containing the blue dot)").shift(3 * DOWN)
        self.play(
            FadeOut(cmartxt, shift = LEFT),
            tarcs[0].animate.shift(-mv),
            impdot.animate.shift(-mv),
        )
        self.wait()
        self.play(Write(expec))
        self.wait()
        hinttxt = Tex("Hint: It's \\\\ not $2\\pi/3$").set_color(YELLOW)
        hinttxt.rotate_about_origin(45 * DEGREES)
        hinttxt.to_corner(UL)
        self.play(Write(hinttxt))
        self.wait(3)


class MotivatingProblem(Scene):
    def construct(self):
        titletxt = Tex(r"An interesting Problem")
        titleline = Underline(titletxt)
        titlegrp = VGroup(titletxt, titleline)
        self.play(Write(titlegrp))
        self.play(titlegrp.animate.to_edge(UP))
        self.wait()
        tri = VMobject(stroke_color = BLUE)
        vera, verb = 6 * RIGHT, 3 * UR + 1.5 * UP
        tri.set_points_as_corners([ORIGIN, vera, verb, ORIGIN])
        tricen = tri.get_center()
        tri.move_to(ORIGIN)
        self.play(Write(tri))
        self.wait()
        dotgrp = VGroup()
        np.random.seed(1)
        for k in range(250):
            x, y = np.random.random(), np.random.random()
            if x + y > 1:
                x, y = 1 - x, 1 - y
            ddot = Dot(x * vera + y * verb, color = GREEN)
            dotgrp.add(ddot)
        dotgrp.shift(-tricen)
        #self.play(LaggedStartMap(FadeIn, dotgrp, run_time = 3))
        self.play(LaggedStart(*[FadeIn(dot, scale = 5) for dot in dotgrp], run_time = 5))
        self.wait()
        verlabs = VGroup()
        verlabs.add(Dot(), Dot(vera), Dot(verb))
        verlabs.add(MathTex("A").next_to(verlabs[0], DOWN, buff = 1 / 8))
        verlabs.add(MathTex("B").next_to(verlabs[1], DOWN, buff = 1 / 8))
        verlabs.add(MathTex("C").next_to(verlabs[2], RIGHT, buff = 1 / 8))
        verlabs.shift(-tricen)
        self.play(Write(verlabs))
        self.wait()
        self.play(*[FadeOut(obj) for obj in dotgrp[1:]])
        self.wait()
        dline = Line(verlabs[0].get_center(), dotgrp[0].get_center(), color = GREEN_E)
        self.play(Create(dline))
        self.wait()
        coord = VGroup(Arrow(verlabs[0].get_center(), verlabs[0].get_center() + RIGHT), Arrow(verlabs[0].get_center(), verlabs[0].get_center() + UP))
        self.play(Write(coord))
        self.wait()
        ml = 3.5 * LEFT
        mr = 2 * RIGHT
        self.play(
            VGroup(verlabs, dline, coord, tri, dotgrp[0]).animate.shift(ml),
            FadeOut(titlegrp, shift = UP)
        )
        dotgrp[1:].shift(ml)
        self.wait()
        verlabs.add(MathTex("P").next_to(dotgrp[0].get_center(), UP))
        integ = MathTex("\\mathbb{E}(AP)", "=", "\\int\\limits_{P \\in \\triangle ABC}", "AP", "\\,\\frac{dydx}{\\triangle ABC}").scale(0.875)
        integ.shift(mr + 2 * UP)
        integ3 = MathTex("\\approx", "\\lim_{n \\to \\infty}", "\\frac{AP_1 + AP_2 + \\cdots + AP_n}{n}").scale(0.875)
        integ3.move_to(integ[1].get_center() + integ3.get_center() - integ3[0].get_center() + 2 * DOWN)
        integ2 = MathTex("=", "\\int\\limits_{P \\in \\triangle ABC}", "\\sqrt{x^2+y^2}", "\\,\\frac{dydx}{\\triangle ABC}").scale(0.875)
        integ2.move_to(integ[1].get_center() + integ2.get_center() - integ2[0].get_center() + 3.5 * DOWN)
        pco = MathTex("(x,y)").scale(0.875).next_to(verlabs[-1], RIGHT)
        self.play(Write(integ), Write(verlabs[-1]))
        verlabs[-1].add_updater(lambda m: m.next_to(dotgrp[0].get_center(), UP))
        dline.add_updater(lambda m: m.become(Line(verlabs[0].get_center(), dotgrp[0].get_center(), color = GREEN_E)))
        self.add(verlabs[-1], dline)
        self.wait()
        animlist = []
        for k in [21, 45, 89, 134, 237, 12, 65, 189, 234, 38, 29, 50]:
            animlist += [Transform(dotgrp[0], dotgrp[k])]
        self.play(
            AnimationGroup(
                Succession(*animlist, lag_ratio = 1.1),
                Write(integ3),
                lag_ratio = 0.5,
                run_time = 8
            )
        )
        self.wait()
        self.play(Write(integ2))
        self.wait()
        self.play(
            FadeOut(Group(verlabs, dline, coord, tri, dotgrp[0]), shift = LEFT),
            FadeOut(integ, shift = UP),
            FadeOut(integ3[0]),
            FadeOut(integ3[1:], shift = RIGHT),
            integ2.animate.move_to(ORIGIN)
        )
        verlabs[-1].clear_updaters()
        dline.clear_updaters()
        self.wait(5)


class IncaseOfASquare(Scene):
    def construct(self):
        integ = VGroup(
            MathTex("\\int\\limits_{P \\in \\triangle ABC}\\text{}"),
            MathTex("\\sqrt{x^2 + y^2}"),
            MathTex("\\,\\frac{dydx}{\\triangle ABC}")
        )
        for txt in integ:
            txt.scale_in_place(0.875)
        integ.arrange(RIGHT)
        self.add(integ)
        self.play(integ.animate.shift(6 * UP))
        self.wait()
        rad, ml = 6, 3.5 * LEFT
        sq = Square(color = BLUE, side_length = rad).shift(ml)
        dotgrp = VGroup()
        a, b = 15, 15
        np.random.seed(1)
        for k in range(a * b):
            xx, yy = 2 * np.random.random() - 1, 2 * np.random.random() - 1
            xx *= rad / 2
            yy *= rad / 2
            dotgrp.add(Dot(xx * RIGHT + yy * UP + ml, color = GREEN))
            dotgrp[-1].save_state()
        dotgrp.arrange_in_grid(a, b).move_to(ORIGIN).shift(-ml)
        np.random.shuffle(dotgrp)
        self.play(Create(sq), LaggedStartMap(FadeIn, dotgrp, run_time = 3))
        self.wait()
        self.play(LaggedStartMap(Restore, dotgrp, path_arc = 75 * DEGREES, rate_func = smooth, run_time = 3))
        self.wait()
        sqintegs = VGroup(
            MathTex("\\int\\limits_{p \\in \\square ABCD}AP\\,\\frac{dydx}{\\square ABCD}"),
            MathTex("=", "\\int\\limits_0^1\\int\\limits_0^1 \\sqrt{x^2 + y^2} \\,dydx"),
            MathTex("=", "\\frac{1}{3}(\\sqrt{2} + \\sinh^{-1}(1))")
        )
        for txt in sqintegs:
            txt.scale_in_place(0.875)
        sqintegs.arrange(DOWN, buff = 1)
        for txt in sqintegs[1:]:
            txt.align_to(sqintegs[0], LEFT)
        sqintegs.shift(-ml)
        labelgrp = VGroup(
            MathTex("A").next_to(sq.get_corner(DL), DOWN),
            MathTex("B").next_to(sq.get_corner(DR), DOWN),
            MathTex("C").next_to(sq.get_corner(UR), UP),
            MathTex("D").next_to(sq.get_corner(UL), UP),
            MathTex("P(x,y)").next_to(dotgrp[0].get_center(), UP),
            Line(sq.get_corner(DL), dotgrp[0].get_center(), color = GREEN)
        )
        txtanimgrp = []
        for txt in sqintegs:
            txtanimgrp += [Write(txt)]
        self.play(
            FadeOut(dotgrp[1:]),
            FadeIn(labelgrp),
            AnimationGroup(*txtanimgrp, lag_ratio = 1.2)
        )
        self.wait()
        self.play(
            FadeOut(VGroup(labelgrp, dotgrp[0], sq), shift = LEFT),
            FadeOut(sqintegs, shift = RIGHT),
            integ.animate.move_to(2 * UP),
        )
        self.wait()

        bub = SVGMobject("/Users/muthuveerappanramalingam/Downloads/ManimInstall/manim_ce/GeometricProbability/assets/Bubbles_thought")
        bub.set_color(WHITE).set_style(fill_color = BLACK, fill_opacity = 3 / 4, stroke_color = WHITE, stroke_width = 3).scale_in_place(2).stretch_in_place(2, 0)
        bub.flip()
        beaidea = Tex("A beautiful \\\\ geometric idea!").set_color(BLUE)
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
        beagrp = add_and_resize(bub.copy(), beaidea)
        beagrp.to_corner(DR)
        bubf = bub.flip()
        solgeom = Tex("Solves a variety \\\\ of problems...").set_color(GREEN)
        solgeomgrp = add_and_resize(bubf, solgeom)
        solgeomgrp.to_corner(DL)
        self.play(Write(beagrp))
        self.wait()
        self.play(Write(solgeomgrp))
        self.wait()
        self.play(
            FadeOut(integ, shift = UP),
            FadeOut(beagrp, shift = RIGHT),
            FadeOut(solgeomgrp, shift = LEFT)
        )
        self.wait()

        tri = VMobject(stroke_color = BLUE)
        vera, verb = 6 * RIGHT, 3 * UR + 1.5 * UP
        tri.set_points_as_corners([ORIGIN, vera, verb, ORIGIN])
        tricen = tri.get_center()
        tri.move_to(ORIGIN)
        self.play(Create(tri))
        dotgrp = VGroup()
        np.random.seed(1)
        for k in range(200):
            x, y = np.random.random(), np.random.random()
            if x + y > 1:
                x, y = 1 - x, 1 - y
            ddot = Dot(x * vera + y * verb, color = GREEN).shift(-tricen)
            dotgrp.add(ddot)
            dotgrp[-1].save_state()
        i = 0
        for dot in dotgrp:
            dot.move_to(tricen)
            dot.shift(5 * RIGHT)
            dot.rotate(about_point = tricen, angle = i)
            i += 5 * DEGREES
        dotgrp.shift(-tricen)
        '''dotgrp[:90].arrange(RIGHT).shift(to_edge(RIGHT))
        dotgrp[91:180].arrange(LEFT).shift(3 * LEFT)
        dotgrp[181:270].arrange(UP).shift(3 * UP)
        dotgrp[271:360].arrange(DOWN).shift(3 * DOWN)'''
        #dotgrp.arrange_in_grid(4, 90).shift(6 * UP)
        self.play(
            LaggedStart(
                #LaggedStartMap(Write, dotgrp),
                #Write(dotgrp),
                LaggedStart(*[FadeIn(dot, scale = 5) for dot in dotgrp]),
                LaggedStartMap(Restore, dotgrp, path_arc = 15 * DEGREES, run_time = 10),
                #lag_ratio = 0.5
            )
        )
        self.wait(3)


class Reduction(Scene):
    def construct(self):
        tri = VMobject(stroke_color = BLUE)
        vera, verb = 6 * RIGHT, 3 * UR + 1.5 * UP
        tri.set_points_as_corners([ORIGIN, vera, verb, ORIGIN])
        tricen = tri.get_center()
        tri.move_to(ORIGIN)
        tric = tri.copy()
        adot, bdot, cdot = Dot(tri.get_start()), Dot(vera - tricen), Dot(verb - tricen)
        randot = Dot()
        labelgrp = VGroup(
            MathTex("A").next_to(adot, DOWN),
            MathTex("B").next_to(bdot, DOWN),
            MathTex("C").next_to(cdot, UP),
        )
        self.play(
            Write(adot), Write(bdot), Write(cdot),
            Write(labelgrp),
            Create(tri)
        )
        self.wait()
        r = ValueTracker(0.5)
        rdeci = DecimalNumber(r.get_value(), num_decimal_places = 2)
        rdeci.add_updater(lambda m: m.set_value(r.get_value()))
        rgrp = VGroup(MathTex("r="), rdeci)
        scarlab = MathTex("\\triangle AB'_rC'_r=r\\cdot \\triangle ABC")
        scarlab.move_to(4 * RIGHT + 2 * UP)
        rgrp[-1].set_color(YELLOW)
        rgrp.arrange(RIGHT)
        rgrp.move_to(4 * RIGHT + 1 * UP)
        #ml = ORIGIN
        def makescaled(t):
            return VGroup(
                tric.copy(),
                bdot.copy().fade(0),
                cdot.copy().fade(0),
            ).scale_about_point(t, adot.get_center())
        def scaleupd(obj):
            obj.become(makescaled(r.get_value()))
        scatri = makescaled(r.get_value())
        bddot, cddot = Dot(scatri[1].get_center()), Dot(scatri[2].get_center())
        bdlab, cdlab = MathTex("B'").next_to(bddot.get_center(), DOWN), MathTex("C'").next_to(cddot.get_center(), UP)
        self.play(
            Write(VGroup(scatri, bddot, cddot, bdlab, cdlab, rgrp, scarlab)),
            tri.animate.fade(0.5)
        )
        bddot.add_updater(lambda m: m.move_to(scatri[1].get_center()))
        cddot.add_updater(lambda m: m.move_to(scatri[2].get_center()))
        bdlab.add_updater(lambda m: m.next_to(bddot.get_center(), DOWN))
        cdlab.add_updater(lambda m: m.next_to(cddot.get_center(), UP))
        self.add(scatri, bddot, cddot, bdlab, cdlab, rgrp)
        self.wait()
        self.play(r.animate.increment_value(0.25), UpdateFromFunc(scatri, scaleupd), run_time = 3)
        self.play(r.animate.increment_value(-0.5), UpdateFromFunc(scatri, scaleupd), run_time = 3)
        self.play(r.animate.increment_value(0.65), UpdateFromFunc(scatri, scaleupd), run_time = 3)
        self.wait()

        pu, pv = 0.3, 0.5
        if pu + pv > 1:
            pu, pv = 1 - pu, 1 - pv
        pdot = Dot(pu * vera + pv * verb).set_color(GREEN).shift(-tricen)
        plab = MathTex("P").next_to(pdot.get_center(), UP)
        self.play(Write(pdot), Write(plab))
        self.wait()
        ml = 3 * LEFT
        self.play(
            VGroup(tri, labelgrp, adot, bdot, cdot, scatri, pdot, plab).animate.shift(ml)
        )
        tric.shift(ml)
        prbexp = VGroup(
            MathTex("\\mathbb{P}(P\\in \\triangle AB'_rC'_r)"),
            MathTex("=\\frac{\\triangle AB'_rC'_r}{\\triangle ABC}"),
            MathTex("=r^2")
        )
        for txt in prbexp:
            txt.scale(0.875)
        prbexp.arrange(RIGHT)
        prbexp.shift(3 * RIGHT + 0.5 * DOWN)
        self.play(Write(prbexp[:-1]))
        self.wait()
        self.play(Write(prbexp[-1]))
        self.wait()

        self.play(
            FadeOut(scarlab, shift = UP),
            rgrp.animate.next_to(tri, DOWN, buff = 0.5),
            FadeOut(prbexp[1]),
            prbexp[2].animate.next_to(prbexp[0], RIGHT)
        )
        self.wait()
        nprbexp = VGroup(prbexp[0], prbexp[2])
        self.play(nprbexp.animate.shift(2 * UR + LEFT))
        self.wait()
        cumustat = Tex(
            "This represents the probability that \\\\ point $P$ is inside $\\triangle AB'_rC'_r$"
        ).scale(0.875)
        #cumustat.move_to(3 * RIGHT + 2 * DOWN)
        cumustat.next_to(nprbexp, UP)
        cumustat[0][40:46].set_color(YELLOW)
        self.play(Write(cumustat))
        self.wait()

        cumueqn = MathTex("\\underline{\\text{CDF:}} \\quad F(r)=r^2")
        cumueqn.next_to(nprbexp, DOWN, buff = 1)
        self.play(Write(cumueqn))
        self.wait()

        pdfeqn = MathTex("\\underline{\\text{PDF:}} \\quad f(r)=F'(r)=2r")
        pdfeqn.next_to(cumueqn, DOWN, buff = 1)
        pdfeqn.shift(0.75 * RIGHT)
        self.play(Write(pdfeqn))
        self.wait()

        densi = Tex(
            "probability or density that \\\\ point $P$ is on side $B'_rC'_r$"
        ).scale(0.875)
        densi.next_to(pdfeqn, DOWN)
        densi[0][32:34].set_color(YELLOW)
        self.play(Write(densi))
        self.wait()

        self.play(
            FadeOut(cumustat, shift = UP),
            FadeOut(nprbexp, shift = UP),
            FadeOut(cumueqn),
            FadeOut(densi, shift = DOWN),
            FadeOut(pdfeqn, shift = DOWN),
        )
        self.wait()

        temparrow = Arrow(ORIGIN, DOWN)
        flowexp = VGroup(
            Tex("Scaling factor, $r$").scale(0.875),
            temparrow.copy(),
            MathTex("F(r)\\text{: }P\\text{ inside }\\triangle AB'_rC'_r").scale(0.875),
            temparrow.copy(),
            MathTex("f(r)\\text{: }P\\text{ on side }B'_rC'_r").scale(0.875),
        )
        flowexp.arrange(DOWN)
        flowexp.shift(3 * RIGHT)
        self.play(
            LaggedStartMap(Write, flowexp, lag_ratio = 0.5),
            r.animate.set_value(0.8), UpdateFromFunc(scatri, scaleupd), run_time = 5   
        )
        self.wait()

        #mr = 7 * RIGHT
        mr = 15 * RIGHT
        self.play(
            FadeOut(flowexp, shift = RIGHT),
            VGroup(tri, labelgrp, adot, bdot, cdot, scatri, pdot, plab, rgrp).animate.shift(mr)
        )
        tric.shift(mr)
        self.wait()

        avgeqns = VGroup(
            Tex("Avg. length of $AP$ given $P$ is inside $\\triangle ABC$"),
            Tex("="),
            MathTex("\\mathbb{P}(P\\text{ is on side }B_rC_r)"),
            MathTex("\\cross"),
            Tex("Avg. length of $AP$ given $P$ is on side $B'_rC'_r$")
        ).arrange(DOWN, buff = 1)
        self.play(Write(avgeqns))
        self.wait()

        sbrace = VGroup(
            MathTex("\\sum_r"),
            Brace(Line(avgeqns[2].get_top(), avgeqns[4].get_bottom()), direction = LEFT).shift(5 * LEFT),
        )
        sbrace[0].next_to(sbrace[1], LEFT)
        self.play(Write(sbrace))
        self.wait()

        expeqn = MathTex("\\mathbb{E}(AP|P \\in \\triangle ABC)=\\int\\limits_0^1 \\mathbb{E}(AP|P\\in B'_rC'_r) \\cdot f(r)\\,dr").scale(0.875)
        self.play(
            FadeOut(avgeqns),
            FadeOut(sbrace),
            FadeIn(expeqn)
        )
        self.wait()
        ml = 8 * LEFT
        self.play(
            VGroup(tri, labelgrp, adot, bdot, cdot, scatri, pdot, plab, rgrp, bddot, cddot, bdlab, cdlab).animate.move_to(ORIGIN).fade(0.875),
            expeqn.animate.move_to(ORIGIN)
        )
        tric.shift(ml)
        self.wait()

        bddot.clear_updaters()
        cdlab.clear_updaters()
        cddot.clear_updaters()
        bdlab.clear_updaters()
        rdeci.clear_updaters()
        self.wait(3)


class GeomOfProblem(Scene):
    def construct(self):
        tri = VMobject(stroke_color = BLUE)
        vera, verb = 6 * RIGHT, 3 * UR + 1.5 * UP
        tri.set_points_as_corners([ORIGIN, vera, verb, ORIGIN])
        tricen = tri.get_center()
        tri.move_to(ORIGIN)
        tric = tri.copy()
        adot, bdot, cdot = Dot(tri.get_start()), Dot(vera - tricen), Dot(verb - tricen)
        randot = Dot()
        labelgrp = VGroup(
            MathTex("A").next_to(adot, DOWN),
            MathTex("B").next_to(bdot, DOWN),
            MathTex("C").next_to(cdot, UP),
        )
        self.play(
            Write(adot), Write(bdot), Write(cdot),
            Write(labelgrp),
            Create(tri)
        )

        x, y = 0.23, 0.33
        if x + y > 1:
            x, y = 1 - x, 1 - y
        r = x + y
        phanttri = tric.copy().scale_about_point(r, adot.get_center())
        phantline = Line(bdot.get_center(), cdot.get_center()).scale_about_point(r, adot.get_center()).fade(1)
        self.add(phantline)
        t = ValueTracker(0.25)
        pdot = Dot(phantline.point_from_proportion(t.get_value()), color = GREEN)
        def makeothers(r):
            pt = phantline.point_from_proportion(r)
            inter = line_intersection(Line(adot.get_center(), pt).get_start_and_end(), Line(bdot.get_center(), cdot.get_center()).get_start_and_end())
            segap = Line(adot.get_center(), pt, color = ORANGE)
            segpp = Line(pt, inter, color = PURPLE)
            pddot = Dot(inter)
            pdlab = MathTex("P'").next_to(pddot, UP)
            return VGroup(segap, segpp, pddot, pdlab)
        seggrp = makeothers(t.get_value())
        #pdot = Dot(phantline.get_start() + t * (phantline.get_end() - phantline.get_start()), color = GREEN)
        bddot, cddot = Dot(phantline.get_start()), Dot(phantline.get_end())
        plab = MathTex("P").next_to(pdot.get_center(), UP)
        bdlab = MathTex("B'").next_to(bddot.get_center(), DOWN)
        cdlab = MathTex("C'").next_to(cddot.get_center(), UP)

        self.play(Write(VGroup(phanttri, pdot, bddot, cddot, plab, bdlab, cdlab)))
        self.wait()
        self.play(
            Write(seggrp),
            tri.animate.fade(0.5)
        )
        self.wait()

        def updothers(obj):
            obj.become(makeothers(t.get_value()))
        
        pdot.add_updater(lambda m: m.become(Dot(phantline.point_from_proportion(t.get_value()), color = GREEN)))
        plab.add_updater(lambda m: m.next_to(pdot.get_center(), UP))
        bdlab.add_updater(lambda m: m.next_to(bddot.get_center(), DOWN))
        cdlab.add_updater(lambda m: m.next_to(cddot.get_center(), UP))
        self.add(pdot, plab, bdlab, cdlab)
        self.wait()

        self.play(
            t.animate.increment_value(0.5),
            UpdateFromFunc(seggrp, updothers),
            rate_func = there_and_back, run_time = 5)
        self.wait()

        trancpy = VGroup(tric.copy(), Dot(seggrp[2].get_center()), Line(adot.get_center(), seggrp[2].get_center(), color = ORANGE))
        self.play(TransformFromCopy(VGroup(phanttri, pdot, seggrp[0]), trancpy, path_arc = 90 * DEGREES), run_time = 3)

        pdot.clear_updaters()
        plab.clear_updaters()
        bdlab.clear_updaters()
        cdlab.clear_updaters()

        scaldist = VGroup(
            MathTex("d(rU, rV)=r \\cdot d(U,V)"),
            Tex("for any two points $U$ and $V$")
        )
        col1, col2 = YELLOW, RED
        print(col1, col2)
        scaldist[0][0][3].set_color(col1)
        scaldist[0][0][6].set_color(col2)
        scaldist[0][0][13].set_color(col1)
        scaldist[0][0][15].set_color(col2)
        scaldist[1][0][15].set_color(col1)
        scaldist[1][0][19].set_color(col2)
        for txt in scaldist:
            txt.scale(0.875)
        scaldist.arrange(DOWN)
        scaldist.move_to(3 * RIGHT - 1.0 * UP)

        gnr = VGroup(
            Tex("Given $r$,"),
            MathTex("\\mathbb{E}(AP|P \\in B'_rC'_r)=r \\cdot \\mathbb{E}(AP'|P'\\in BC)")
        )
        for txt in gnr:
            txt.scale(0.875)
        gnr.arrange(DOWN)
        gnr.move_to(3 * RIGHT - 1.0 * DOWN)

        self.play(
            Write(scaldist),
            VGroup(phanttri, pdot, bddot, cddot, plab, bdlab, cdlab, seggrp, labelgrp, tri, adot, bdot, cdot, trancpy).animate.shift(3 * LEFT)
        )
        self.wait()
        self.play(Write(gnr))
        self.wait()

        self.play(
            FadeOut(VGroup(phanttri, pdot, bddot, cddot, plab, bdlab, cdlab, seggrp, labelgrp, tri, adot, bdot, cdot, trancpy), shift = LEFT),
            FadeOut(scaldist),
            FadeOut(gnr[0]),
            gnr[1].animate.move_to(ORIGIN)
        )
        self.wait()

        firstrel = MathTex("\\mathbb{E}(AP|P\\in \\triangle ABC)=\\int\\limits_0^1 \\mathbb{E}(AP|P \\in B'_rC'r)\\cdot f(r) \\, dr")
        firstrel.scale(0.875)
        firstrel.shift(2 * UP)
        remi = MathTex("f(r)=2r").scale(0.75)
        remirec = SurroundingRectangle(remi)
        remigrp = VGroup(remi, remirec).to_corner(UR)

        self.play(
            Write(firstrel),
            Write(remigrp)
        )
        self.wait()

        finalrel = MathTex("\\mathbb{E}(AP|P\\in \\triangle ABC)=\\int\\limits_0^1 r\\cdot \\mathbb{E}(AP'|P' \\in BC)\\cdot 2r \\, dr")
        finalrel.scale(0.875)
        finalrel.shift(2 * DOWN)

        self.play(Write(finalrel))
        self.wait()
        self.play(
            LaggedStart(
                FadeOut(firstrel, shift = UP),
                FadeOut(gnr[1], shift = UP),
                finalrel.animate.move_to(ORIGIN)
            )
        )
        self.wait()

        firstsimp = MathTex("\\mathbb{E}(AP|P\\in \\triangle ABC)=\\int\\limits_0^1 \\mathbb{E}(AP'|P' \\in BC)\\cdot 2r^2 \\, dr").scale(0.875)
        self.play(
            FadeOut(remigrp),
            TransformMatchingTex(finalrel, firstsimp)
        )
        self.wait()

        indpr = Tex("independent of $r$!!").shift(2 * UR)
        arr = Arrow(firstsimp[0][20].get_top(), indpr.get_bottom(), color = BLUE)

        self.play(Write(indpr), Create(arr))
        self.wait()

        secondsimp = MathTex("\\mathbb{E}(AP|P\\in \\triangle ABC)=\\mathbb{E}(AP'|P' \\in BC)\\int\\limits_0^1 2r^2 \\, dr").scale(0.875)
        self.play(
            TransformMatchingTex(firstsimp, secondsimp),
            FadeOut(VGroup(indpr, arr))
        )
        self.wait()

        finalsimp = MathTex("\\mathbb{E}(AP|P\\in \\triangle ABC)", "=", "\\frac{2}{3}\\cdot", "\\mathbb{E}(AP'|P' \\in BC)").scale(0.875)
        b1 = BraceText(finalsimp[0], "2-D integral")
        b2 = BraceText(finalsimp[-1], "1-D integral")
        self.play(
            FadeOut(secondsimp, shift = UP),
            FadeIn(finalsimp, shift = UP),
        )
        self.wait()
        self.play(Write(b1))
        self.wait()
        self.play(TransformFromCopy(b1, b2, path_arc = 120 * DEGREES), run_time = 3)
        self.wait()
        self.play(
            FadeOut(VGroup(b1, b2), shift = UP),
            finalsimp.animate.shift(2 * UP),
            FadeOut(finalsimp[:-1], shift = UP),
            finalsimp[-1].animate.to_edge(LEFT).shift(2 * UP)
        )
        self.wait()
        othergp = VGroup(
            Tex("perimater, $p=a+b+c$"),
            Tex("semi-perimeter, $s=(a+b+c)/2$"),
            Tex("length of altitude to side $BC$, $h_a$")
        )
        for txt in othergp:
            txt.scale(0.875)
        othergp.arrange(DOWN)
        othergp.shift(0.5 * DOWN)
        expec = MathTex("=\\frac{b+c}{4}\\left(1+\\left(\\frac{b-c}{a}\\right)^2\\right)+\\frac{h_a^2}{2}\\csch^{-1}\\left(\\frac{p}{a}\\cdot\\frac{s-a}{p-a}\\right)").scale(0.875)
        expec.next_to(finalsimp[-1], RIGHT)
        self.play(
            AnimationGroup(
                Write(expec),
                Write(othergp),
                lag_ratio = 0.75
            )
        )
        self.wait(3)


class ThisToThat(Scene):
    def construct(self):
        tri = VMobject(stroke_color = BLUE)
        vera, verb = 6 * RIGHT, 3 * UR + 1.5 * UP
        tri.set_points_as_corners([ORIGIN, vera, verb, ORIGIN])
        tricen = tri.get_center()
        tri.move_to(ORIGIN)
        tric = tri.copy()

        areagp, sidegp = VGroup(), VGroup()
        bc = Line(vera, verb)
        c1, c2 = YELLOW, ORANGE
        mr = 3.5 * RIGHT
        for k in range(100):
            x, y, z = np.random.random(), np.random.random(), np.random.random()
            if x + y > 1:
                x, y = 1 - x, 1 - y
            arpos, sipos = x * vera + y * verb, bc.point_from_proportion(z)
            ardot, sidot = Dot(arpos, color = c1), Dot(sipos, color = c2)
            ardot.shift(-tricen - mr)
            sidot.shift(-tricen + mr)
            areagp.add(ardot)
            sidegp.add(sidot)
            areagp[-1].save_state()
            sidegp[-1].save_state()
        tri.shift(mr)
        tric.shift(-mr)
        for obj in areagp:
            obj.move_to(ORIGIN)
        for obj in sidegp:
            obj.move_to(ORIGIN)
        dotgrp = VGroup(Dot(), Dot(vera), Dot(verb))
        labgrp1 = VGroup(
            MathTex("A").next_to(dotgrp[0], DOWN),
            MathTex("B").next_to(dotgrp[1], DOWN),
            MathTex("C").next_to(dotgrp[2], RIGHT),
        ).shift(-tricen)
        labgrp2 = labgrp1.copy()
        labgrp1.shift(mr)
        labgrp2.shift(-mr)
        self.play(
            Create(tri), Create(tric),
            FadeIn(areagp), FadeIn(sidegp),
            Write(labgrp1), Write(labgrp2)
        )
        self.wait()
        txtgrp = MathTex("\\mathbb{E}(AP|P\\in \\triangle ABC)", "=", "\\frac{2}{3}", "\\cross", "\\mathbb{E}(AP|P \\in BC)")
        txtgrp.to_edge(UP)
        self.play(
            LaggedStart(*[Restore(obj, path_arc = 60 * DEGREES) for obj in areagp]),
            Write(txtgrp),
            LaggedStart(*[Restore(obj, path_arc = 60 * DEGREES) for obj in sidegp]),
        )
        self.wait()
        self.wait()


class OtherProblems(Scene):
    def construct(self):
        rad = 3
        sq, cir = Square(side_length = 2 * rad, color = BLUE), Circle(radius = rad, color = BLUE)
        sqgp, cirgp = VGroup(), VGroup()
        mr = (rad + 1 / 2) * RIGHT
        #col1, col2 = GREY, GREY
        col1, col2 = GREEN, GREEN
        for k in range(100):
            u, v, x, y = np.random.random(), np.random.random(), np.random.random(), np.random.random()
            u = rad * (2 * u - 1)
            v = rad * (2 * v - 1)
            x = rad * (2 * x - 1)
            y = rad * (2 * y - 1)
            fseg = Line(x * RIGHT + y * UP, u * RIGHT + v * UP, color = col1)
            fseg.shift(mr)
            sqgp.add(fseg)
            sqgp[-1].save_state()
            u, v, x, y = (u / rad + 1) / 2, (v / rad + 1) / 2, (x / rad + 1) / 2, (y / rad + 1) / 2
            dot1 = Dot(rad * np.sqrt(u) * RIGHT).rotate_about_origin(x * 360 * DEGREES)
            dot2 = Dot(rad * np.sqrt(v) * RIGHT).rotate_about_origin(y * 360 * DEGREES)
            radline = Line(dot1.get_center(), dot2.get_center(), color = col2)
            radline.shift(-mr)
            cirgp.add(radline)
            cirgp[-1].save_state()
        sq.shift(mr)
        cir.shift(-mr)
        self.play(
            Write(sq), Write(cir)
        )
        self.wait()
        for obj in cirgp:
            obj.move_to(ORIGIN)
        for obj in sqgp:
            obj.move_to(ORIGIN)
        #sqgp.shift(10 * LEFT + 5 * DOWN).fade(1)
        #cirgp.shift(10 * RIGHT + 5 * UP).fade(1)
        sqgp.to_corner(UR).fade(1)
        cirgp.to_corner(DL).fade(1)
        self.play(
            LaggedStart(*[Restore(obj, path_arc = 180 * DEGREES) for obj in cirgp]),
            LaggedStart(*[Restore(obj, path_arc = 180 * DEGREES) for obj in sqgp]),
        )
        self.wait()
        expexp = VGroup(
            Tex("$\\mathbb{E}($random segment in a circle$)$?").scale(0.875),
            Tex("$\\mathbb{E}($random segment in a square$)$?").scale(0.875)
        )
        expexp[0].next_to(cir, DOWN)
        expexp[1].next_to(sq, DOWN)
        self.play(Write(expexp))
        self.wait(3)


if __name__ == "__main__":
    module_name = os.path.abspath(__file__)
    #output_location = "C:\ManimCE\media"
    #clear_cmd = "cls"
    #command_A = "manim " + module_name + " " + "RealCase" + " " + "-pql -n 42" + " --media_dir " + output_location
    #command_A = "manim " + module_name + " --media_dir " + output_location + " -pqh"
    command_A = "manim "+ "-pql" + " " + module_name + " " + "InterviewQuestion" + " -n 0,5"
    #command_A = "manim "+ "-pqh" + " " + module_name
    #os.system(clear_cmd)
    os.system(command_A)
