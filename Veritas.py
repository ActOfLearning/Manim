from manim import *
import os
from manim.mobject.svg.style_utils import SVG_DEFAULT_ATTRIBUTES

from manim.utils.file_ops import add_version_before_extension
from numpy import dot


class TuesdayGirl(Scene):
    def construct(self):
        strwid = 1 / 2
        bb = Square(fill_opacity = 1, fill_color = BLUE, stroke_width = strwid)
        bg = Square(fill_opacity = 1, fill_color = color_gradient([BLUE, BLUE, GREEN, GREEN], 10), stroke_width = strwid)
        gb = Square(fill_opacity = 1, fill_color = color_gradient([GREEN, GREEN, BLUE, BLUE], 10), stroke_width = strwid)
        gg = Square(fill_opacity = 1, fill_color = GREEN, stroke_width = strwid)
        childgrp = VGroup(bb, bg, gb, gg)
        childgrp.arrange_in_grid(2, 2, buff = 1 / 8)
        daychildgrp = VGroup()
        for i in range(49):
            daychildgrp.add(childgrp.copy())
        daychildgrp.arrange_in_grid(7, 7, buff = 1 / 4)
        daychildgrp.shift(0.25 * UP)

        tuesgrp, nontuesgrp = VGroup(), VGroup()
        tuesboysgrp = VGroup()
        errgrp1, errgrp2 = VGroup(), VGroup()
        for i in range(49):
            if i % 7 == 2 or i // 7 == 4:
                tuesgrp.add(daychildgrp[i])
                tuesboysgrp.add(daychildgrp[i][0])
                if i % 7 == 2 and not(i // 7 == 4):
                    errgrp1.add(daychildgrp[i][1])
                if not(i % 7 == 2) and i // 7 == 4:
                    errgrp2.add(daychildgrp[i][2])
            else:
                nontuesgrp.add(daychildgrp[i])
        
        veritas = Tex("Veritasium").scale(2)
        veriline = Underline(veritas)
        verigrp = VGroup(veritas, veriline)
        self.play(Write(verigrp))
        self.wait()

        tuestitle = Tex("Tuesday Girl Problem").scale(1.25)
        tuestitle[0][:7].set_color(YELLOW)
        tuestitle[0][7:11].set_color(GREEN)
        titleline = Underline(tuestitle)
        titlegrp = VGroup(tuestitle, titleline)
        self.play(ReplacementTransform(verigrp, titlegrp))
        self.wait()
        self.play(titlegrp.animate.to_edge(UP), buff = 1 / 16)
        self.wait()

        consitxt = Tex("Considering a family of two children,\\\\ with ", "boys", " and ", "girls", " equally likely")
        consitxt[1].set_color(BLUE)
        consitxt[3].set_color(GREEN)
        consitxt.shift(1.5 * UP)
        consitxt2 = Tex("child's gender is ", "independent", " of the day of \\\\ week on which it is born")
        consitxt2[1].set_color(RED)
        self.play(Write(consitxt))
        self.play(Write(consitxt2))
        #self.play(Write(consitxt), Write(consitxt2))
        self.wait()

        legen = childgrp.copy()
        legen.scale_in_place(5 / 8)
        legen = VGroup(legen[0], legen[2], legen[1], legen[3])
        a, b, c = 2, 2, 1.5
        legen[0].move_to(-a * RIGHT - b * UP)
        legen[1].move_to(a * RIGHT - b * UP)
        legen[2].move_to(-a * RIGHT - (b + c) * UP)
        legen[3].move_to(a * RIGHT - (b + c) * UP)
        leglabels = VGroup(
            Tex("boy", "boy"),
            Tex("boy", "girl"),
            Tex("girl", "boy"),
            Tex("girl", "girl")
        )
        for obj in leglabels:
            obj.arrange(DOWN, buff = 1 / 16)
        for j in range(4):
            col1 = BLUE if j % 2 == 0 else GREEN
            col2 = BLUE if j // 2 == 0 else GREEN
            leglabels[j][0].set_color(col2)
            leglabels[j][1].set_color(col1)
        for j in range(4):
            leglabels[j].next_to(legen[j], LEFT)
        self.play(Write(leglabels), Create(legen))
        self.wait()

        self.play(
            FadeOut(VGroup(consitxt, consitxt2, legen, leglabels)),
        )
        self.wait()

        prob1 = MathTex("\\mathbb{P}(\\text{both girls|atleast one girl})")
        prob1[0][6:11].set_color(GREEN)
        prob1[0][22:26].set_color(GREEN)
        chil = childgrp.copy().scale(0.875)
        self.play(Write(prob1))
        self.wait()

        chil.to_edge(DR)
        self.play(
            prob1.animate.shift(2 * DOWN + 3 * LEFT),
            LaggedStartMap(GrowFromCenter, chil)
        )
        self.wait()

        res1 = MathTex("=\\frac{1}{3}").next_to(prob1, RIGHT)
        self.play(
            FadeOut(chil[0]),
            Write(res1)
        )
        self.wait()

        prob2 = MathTex("\\mathbb{P}(\\text{both girls|elder child is a girl})").shift(1.5 * UP)
        prob2[0][6:11].set_color(GREEN)
        prob2[0][25:29].set_color(GREEN)
        #prob2[0][28:35].set_color(YELLOW)
        self.play(
            FadeOut(chil[1:], shift = RIGHT),
            FadeOut(res1, shift = RIGHT),
            FadeOut(consitxt),
            ReplacementTransform(prob1, prob2)
        )
        self.wait()

        chil2 = childgrp.copy().shift(1.5 * DOWN)
        self.play(
            LaggedStart(*[FadeIn(obj, shift = UL) for obj in chil2])
        )
        oldyo = Tex("Elder", "Younger")
        stc = chil2[1].get_center()
        stu = chil2[1].get_corner(UL)
        stl = chil2[1].get_corner(DR)
        stu = (stc + stu) / 2
        stl =(stc + stl) / 2
        twoarr = VGroup(
            Arrow(stu, stu + 2.5 * RIGHT),
            Arrow(stl, stl + 1.5 * RIGHT)
        )
        oldyo[0].next_to(twoarr[0], RIGHT)
        oldyo[1].next_to(twoarr[1], RIGHT)
        self.play(Create(twoarr), Write(oldyo))
        self.wait()

        self.play(
            FadeOut(oldyo, shift = RIGHT),
            FadeOut(twoarr, shift = RIGHT),
        )
        self.wait()
        res2 = MathTex("=\\frac{1}{2}").next_to(prob2)
        self.play(FadeOut(chil2[0]), FadeOut(chil2[2]), Write(res2))
        self.wait()

        prob4 = MathTex("\\mathbb{P}(\\text{both girls|elder girl born on a tuesday})").to_edge(UP)
        prob4[0][6:11].set_color(GREEN)
        prob4[0][17:21].set_color(GREEN)
        prob4[0][28:35].set_color(YELLOW)
        self.play(
            FadeOut(titlegrp, shift = UP),
            FadeOut(chil2[1]), FadeOut(chil2[3]),
            FadeOut(res2, shift = UP),
            ReplacementTransform(prob2, prob4)
        )
        self.wait()

        weekdays = Tex("Su", "Mo", "Tu", "We", "Th", "Fr", "Sa").set_color_by_gradient(YELLOW, RED, PINK, PURPLE)
        wkdys = weekdays.copy().scale(0.75)
        wkdysc = wkdys.copy().rotate_in_place(90 * DEGREES)
        daymat = daychildgrp.copy().scale(0.171875)
        daymat.shift(0.5 * DOWN)
        for i in range(7):
            wkdys[i].next_to(daymat[42 + i], DOWN)
            wkdysc[i].next_to(daymat[42 - 7 * i], LEFT)
        eldtues, neldtues = VGroup(), VGroup()
        daylist = [VGroup() for i in range(7)]
        for i in range(49):
            if i // 7 == 4:
                eldtues.add(VGroup(daymat[i][1], daymat[i][3]))
                #neldtues.add(VGroup(daymat[i][0], daymat[i][2]))
                daylist[4].add(VGroup(daymat[i][0], daymat[i][2]))
            else:
                daylist[i // 7].add(daymat[i])
        wktxt = VGroup(
            Tex("Elder child").rotate_in_place(90 * DEGREES),
            Tex("Younger child")
        )
        wktxt[1].next_to(wkdys, DOWN)
        wktxt[0].next_to(wkdysc, LEFT)
        '''for obj in daymat:
            obj.save_state()
            xrand, yrand, rot = np.random.random(), np.random.random(), np.random.random()
            #obj.move_to(7 * (2 * xrand - 1) * RIGHT + 4 * (2 * yrand - 1) * UP)
            obj.shift(7 * (2 * xrand - 1) * RIGHT)
            obj.rotate_in_place(rot * 720 * DEGREES)'''
        self.play(LaggedStartMap(GrowFromCenter, daymat), run_time = 2)
        #self.play(LaggedStartMap(Restore, daymat), run_time = 3)
        self.play(Write(wkdys), Write(wkdysc), Write(wktxt))
        self.wait()
        self.play(LaggedStartMap(FadeOut, daylist), run_time = 2)
        self.wait()
        res2.next_to(prob4, RIGHT)
        self.play(FadeIn(res2))
        self.wait()
        '''self.play(LaggedStartMap(FadeOut, neldtues), run_time = 1.5)
        self.wait()'''

        prob3 = MathTex("\\mathbb{P}(\\text{both girls|atleast one girl on tuesday})").to_edge(UP)
        prob3[0][6:11].set_color(GREEN)
        prob3[0][12:19].set_color(RED)
        prob3[0][22:26].set_color(GREEN)
        prob3[0][28:35].set_color(YELLOW)
        self.play(
            FadeOut(VGroup(wktxt[0], wkdysc), shift = LEFT),
            FadeOut(VGroup(wktxt[1], wkdys), shift = DOWN),
            FadeOut(VGroup(res2, eldtues)),
            ReplacementTransform(prob4 , prob3)
        )
        self.wait()

        ques1 = Tex("Isn't it similar to \\\\ the previous question?")
        ques2 = Tex("Anyway, weekdays \\\\ doesn't matter")
        bub = SVGMobject("/Users/muthuveerappanramalingam/Downloads/ManimInstall/manim_ce/GeometricProbability/assets/Bubbles_thought")
        bub.set_color(WHITE).set_style(fill_color = BLACK, fill_opacity = 3 / 4, stroke_color = WHITE, stroke_width = 3).scale_in_place(2).stretch_in_place(2, 0)
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
        ques1grp = add_and_resize(bub.copy().set_color_by_gradient(BLUE, YELLOW), ques1)
        ques2grp = add_and_resize(bub.copy().flip().set_color_by_gradient(YELLOW, RED), ques2)
        ques1grp.to_corner(DR)
        ques2grp.to_corner(DL)
        self.play(Write(ques1grp))
        self.play(Write(ques2grp))
        self.wait()

        daychildgrp.scale(1 / 4 - 1 / 64).shift(2 * LEFT)
        weekdays.scale(0.875)
        weekdaysc = weekdays.copy().rotate_about_origin(90 * DEGREES)
        for i in range(7):
            weekdays[i].next_to(daychildgrp[42 + i], DOWN)
            weekdaysc[i].next_to(daychildgrp[42 - 7 * i], LEFT)
        for obj in daychildgrp:
            obj.save_state()
            sc, rot = np.random.random(), np.random.random()
            obj.scale_in_place(sc)
            obj.rotate_in_place(rot * rot * 360 * DEGREES)
            obj.shift(2 * RIGHT)
        self.play(
            FadeOut(prob3, shift = UP),
            FadeOut(ques1grp, shift = RIGHT),
            FadeOut(ques2grp, shift = LEFT),
            LaggedStart(*[FadeIn(obj, scale = 1 + np.random.random()) for obj in daychildgrp], run_time = 1.5)
            #LaggedStartMap(GrowFromCenter, daychildgrp, run_time = 1.5)
        )
        #self.play(LaggedStartMap(Restore, daychildgrp))
        self.play(LaggedStart(*[Restore(obj, path_arc = np.random.random() * 135 * DEGREES) for obj in daychildgrp]), run_time = 3)
        self.play(Write(weekdays), Write(weekdaysc))
        self.wait()

        finaltxt = VGroup(
            Tex("$\\cross$", "both non-tuesdays"),
            Tex("$\\cross$", "both boys"),
            Tex("$\\cross$", "no tuesday girls")
        )
        finaltxt.arrange(DOWN, buff = 1)
        finaltxt[1].align_to(finaltxt[0], LEFT)
        finaltxt[2].align_to(finaltxt[0], LEFT)
        finaltxt.shift(4.25 * RIGHT)
        finaltxt[0][0].set_color(RED)
        finaltxt[1][0].set_color(RED)
        finaltxt[2][0].set_color(RED)
        finaltxt[0][1][8:].set_color(YELLOW)
        finaltxt[1][1][4:].set_color(BLUE)
        finaltxt[2][1][9:].set_color(GREEN)
        finaltxt[2][1][2:9].set_color(YELLOW)
        self.play(
            LaggedStartMap(FadeOut, nontuesgrp),
            Write(finaltxt[0])
        )
        self.wait()
        self.play(
            LaggedStartMap(FadeOut, tuesboysgrp),
            Write(finaltxt[1])
        )
        self.wait()
        self.play(
            LaggedStartMap(FadeOut, errgrp1),
            LaggedStartMap(FadeOut, errgrp2),
            Write(finaltxt[2])
        )
        self.wait()

        self.play(LaggedStart(*[FadeOut(obj, shift = DOWN) for obj in finaltxt]))
        self.wait()
        probtxt = VGroup(
            Tex("$\\mathbb{P}$(both girls|atleast \\\\ one girl on a tuesday)"),
            MathTex("="),
            MathTex("13", "\\over", "27")
        )
        probtxt.arrange(DOWN)
        probtxt.shift(4.25 * RIGHT)
        probtxt[1].next_to(probtxt[2], LEFT)
        probtxt[0][0][6:11].set_color(GREEN)
        probtxt[0][0][12:19].set_color(RED)
        probtxt[0][0][22:26].set_color(GREEN)
        probtxt[0][0][29:36].set_color(YELLOW)
        probtxt[2][0].set_color(GREEN)
        probtxt[2][-1].set_color_by_gradient(BLUE, GREEN, BLUE, GREEN)
        self.play(Write(probtxt[0]))
        self.play(Write(probtxt[1]), Write(probtxt[2]))
        self.wait(5)


if __name__ == "__main__":
    module_name = os.path.abspath(__file__)
    #output_location = "C:\ManimCE\media"
    #clear_cmd = "cls"
    #command_A = "manim " + module_name + " " + "RealCase" + " " + "-pql -n 42" + " --media_dir " + output_location
    #command_A = "manim " + module_name + " --media_dir " + output_location + " -pqh"
    #command_A = "manim "+ "-pql" + " " + module_name + " " + "TuesdayGirl" + " -n 0"
    command_A = "manim "+ "-pqh" + " " + module_name
    #os.system(clear_cmd)
    os.system(command_A)
