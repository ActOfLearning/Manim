from manimlib.imports import *
import numpy as np

class LearningScene(GraphScene, MovingCameraScene):
    CONFIG = {
        "x_min": -10,
        "x_max": 10,
        "x_axis_width": 20,
        "x_tick_frequency": 1,
        "x_leftmost_tick": None,  # Change if different from x_min
        "x_labeled_nums": None,
        "x_axis_label": "",
        "y_min": -6,
        "y_max": 6,
        "y_axis_height": 12,
        "y_tick_frequency": 1,
        "y_bottom_tick": None,  # Change if different from y_min
        "y_labeled_nums": None,
        "y_axis_label": "",
        "axes_color": GREY,
        "graph_origin": ORIGIN,
        "exclude_zero_label": True,
        "default_graph_colors": [BLUE, GREEN, YELLOW],
        "default_derivative_color": GREEN,
        "default_input_color": YELLOW,
        "default_riemann_start_color": BLUE,
        "default_riemann_end_color": GREEN,
        "area_opacity": 0.8,
        "num_rects": 50,
        "dot_kwargs": {
            "radius": 0.05,
        },
        "line_kwargs": {
            "stroke_width": 2,
        },
        "fill_triangle_kwargs": {
            # "fill_color": BLUE,
            "fill_opacity": .5,
            "stroke_width": 0,
        },
    }

    def setup(self):
        GraphScene.setup(self)
        MovingCameraScene.setup(self)
    
    def construct(self):
        self.setup_axes()
        #self.remove(self.axes)
        self.play(FadeOut(self.axes), run_time = 0.01)
        f_x = VGroup(
            TexMobject("\\displaystyle\\int\\limits_{0}^b"),
            TexMobject("f(x)"),
            TexMobject("\\,dx")
            ).arrange(direction = RIGHT, buff = 0.1)
        integ = VGroup(
            TexMobject("\\displaystyle\\int\\limits_{0}^b"),
            TexMobject("\\sqrt{a^2-x^2}").set_color(BLUE),
            TexMobject("\\,dx")
        ).arrange(direction = RIGHT, buff = 0.1)
        self.play(Write(integ), run_time = 5)
        self.wait()
        result = VGroup(
            TexMobject("\\frac{1}{2}b\\sqrt{a^2-b^2}").set_color(YELLOW),
            TexMobject("+"),
            TexMobject("\\frac{1}{2}a^2\\sin^{-1}\\frac{b}{a}").set_color(GREEN)
        ).arrange(direction = RIGHT, buff = 0.125)
        eq_sign = TexMobject("=").move_to(LEFT)
        self.play(
            integ.next_to, eq_sign, {"direction": LEFT, "buff": 0.25})
        result.next_to(eq_sign, direction = RIGHT, buff = 0.125)
        self.wait()
        self.play(Write(eq_sign), Write(result))
        self.wait()
        self.play(*[FadeOut(i) for i in [eq_sign, result]], integ.move_to, ORIGIN)
        self.wait()
        tbubble = ThoughtBubble()
        tbubble.pin_to(integ)
        seensom = TextMobject("... solved this in \\\\ high school...")
        bub_text = seensom.copy()
        tbubble.add_content(bub_text)
        tbubble.resize_to_content()
        bub_grp = VGroup(tbubble, bub_text)
        self.play(Write(bub_grp))
        self.wait()
        trigsub = TextMobject("Trigonometric \\\\ substitution!!")
        tbubble.add_content(trigsub)
        self.play(Transform(bub_text, trigsub))
        self.wait()
        self.play(FadeOut(bub_grp))
        self.wait()
        self.play(integ.shift, 2.5 * UP)
        self.wait()
        solplan = VGroup(
            TextMobject("$\\bullet$ make the substitution $x=a \\sin \\theta$"),
            TextMobject("$\\bullet$ evaluate $dx = a \\cos \\theta d\\theta$"),
            TextMobject("$\\bullet$ change the limits of the integral to $0 \\leq \\theta \\leq \\sin^{-1}(x/a)$"),
            TextMobject("$\\bullet$ solve the resulting integral"),
            TextMobject("$\\bullet$ DONE!!")
        ).arrange(direction = DOWN, buff = 0.5)
        for text in solplan:
            text.scale_in_place(0.75)
            text.align_to(integ, LEFT)
        solplan[-1].scale_in_place(20 / 19)
        solplan.shift(4 * LEFT + DOWN)
        for text in solplan:
            self.play(Write(text))
            self.wait()
        self.play(FadeOut(solplan))
        self.wait()
        uldot, urdot = Dot().to_corner(DL), Dot().to_corner(DR)
        ulbub, urbub = SpeechBubble(), SpeechBubble()
        ulbub.pin_to(uldot)
        urbub.pin_to(urdot)
        ultxt = TextMobject("solve the problem..")
        urtxt = TextMobject("Is there a better way?")
        ulbub.add_content(ultxt)
        urbub.add_content(urtxt)
        ulbub.resize_to_content()
        urbub.resize_to_content()
        ulgrp = VGroup(ulbub, ultxt)
        urgrp = VGroup(urbub, urtxt)
        self.play(Write(ulgrp))
        self.wait()
        self.play(Write(urgrp))
        self.wait()
        newpers = TextMobject("A different persective!!").scale(1.25).set_color(YELLOW).move_to(integ.get_center())
        newpers.shift(3 * DOWN)
        self.play(FadeOut(ulgrp), FadeOut(urgrp), FadeInFromDown(newpers))
        self.wait()
        self.play(FadeOut(newpers))
        self.wait()
        f_x.move_to(integ)
        self.remove(integ)
        self.play(ReplacementTransform(integ.copy(), f_x))
        self.wait()
        self.play(
            self.camera_frame.move_to, 2.5 * UP + 5 * RIGHT,
            f_x.move_to, 5.25 * UP + 5 * RIGHT
        )
        self.play(FadeIn(self.axes))
        self.wait()
        func = self.get_graph(lambda t: (11 * t ** 3 - 139 * t ** 2 + 6 * 83 * t) / 180, x_min = -2, x_max = 9)
        self.play(Write(func))
        self.wait()
        rect_list = self.get_riemann_rectangles_list(
            func, 7, 
            max_dx = 1,
            x_min = 0,
            x_max = 8,
        )
        flat_graph = self.get_graph(lambda t : 0)
        rects = self.get_riemann_rectangles(
            flat_graph, x_min = 0, x_max = 8, dx = 1.0
        )
        for new_rects in rect_list:
            new_rects.set_fill(opacity = 0.8)
            rects.align_submobjects(new_rects)
            for alt_rect in rects[::2]:
                alt_rect.set_fill(opacity = 0)
            self.play(Transform(
                rects, new_rects,
                run_time = 2,
                lag_ratio = 0.5
            ))
        self.play(FadeOut(rects), FadeOut(func))
        self.wait()
        cfunc = self.get_graph(lambda t: np.sqrt(25 - t ** 2), x_min = 0, x_max = 5, color = BLUE)
        self.remove(f_x)
        integ.move_to(f_x.get_center())
        self.play(ReplacementTransform(f_x.copy(), integ))
        yeq = VGroup(
            TexMobject("y=\\sqrt{a^2-x^2}").move_to(2.25 * RIGHT + 2 * UP),
            TexMobject("x^2+y^2=a^2").move_to(2.25 * RIGHT + 2 * UP)
            )
        self.play(FadeOut(func), Write(cfunc), integ.shift, 8.5 * LEFT + 3 * DOWN, Write(yeq[0]), self.camera_frame.shift, 4 * LEFT)
        self.wait()
        self.play(Transform(yeq[0], yeq[1]))
        self.wait()
        self.play(FadeOut(yeq))
        self.wait()
        v_line = DashedLine(3 * RIGHT + 0.5 * DOWN, 3 * RIGHT + 6 * UP).set_color(GRAY)
        limlabels = VGroup(TexMobject("x=0"), TexMobject("x=b"))
        for text in limlabels:
            text.rotate_in_place(PI / 2).scale_in_place(0.75)
        limlabels[1].next_to(v_line.get_end(), direction = DL, buff = 0.1)
        limlabels[0].move_to(limlabels[1].get_center())
        limlabels[0].shift(3 * LEFT)
        limlabels.add(v_line)
        self.play(Write(limlabels))
        self.wait()
        secarea = VMobject(fill_color = GREEN, fill_opacity = 0.75, stroke_width = 0)
        secarea_points = [self.input_to_graph_point(i, cfunc) for i in np.arange(0, v_line.get_start()[0] + 0.1, 0.1)]
        segarea = VMobject(fill_color = BLUE, fill_opacity = 0.75, stroke_width = 0)
        segarea.set_points_as_corners(secarea_points + [3 * RIGHT, ORIGIN])
        secarea_points.append(ORIGIN)
        secarea.set_points_as_corners(secarea_points)
        triarea = VMobject(fill_color = YELLOW, fill_opacity = 0.75, stroke_width = 0)
        triarea.set_points_as_corners([ORIGIN, 3 * RIGHT, 3 * RIGHT + 4 * UP])
        area_gp = VGroup(secarea, triarea)
        rad_brace = Brace(Line(ORIGIN, 5 * UP), direction = LEFT)
        rad_label = TexMobject("a").scale(0.75).next_to(rad_brace, direction = LEFT, buff = 0.1)
        self.play(Write(segarea), Write(rad_brace), Write(rad_label))
        self.wait()
        self.remove(segarea)
        self.play(ReplacementTransform(segarea.copy(), area_gp))
        self.wait()
        b_brace = Brace(Line(ORIGIN, 3 * RIGHT), direction = DOWN)
        sq_brace = Brace(Line(3 * RIGHT, 3 * RIGHT + 4 * UP), direction = RIGHT)
        b_label = TexMobject("b").scale(0.75).next_to(b_brace, direction = DOWN, buff = 0.1)
        sq_label = TexMobject("\\sqrt{a^2-b^2}").scale(0.75).rotate(-PI / 2).next_to(sq_brace, direction = RIGHT, buff = 0.1)
        b_line = DashedLine(4 * UP, 4 * UP + 3 * RIGHT).set_color(GREEN)
        b_sec = b_label.copy().next_to(b_line, direction = UP, buff = 0.1)
        self.play(*[Write(i) for i in [b_brace, sq_brace, b_label, sq_label, b_line, b_sec]])
        self.wait()
        tri_gp = VGroup(triarea, b_brace, b_label, sq_brace, sq_label)
        sec_gp = VGroup(secarea, b_line, b_sec, rad_brace, rad_label)
        feq_sign = eq_sign.copy().next_to(integ, direction = RIGHT, buff = 0.25)
        self.play(*[FadeOut(i) for i in [self.axes, limlabels, cfunc]], tri_gp.shift, 3 * RIGHT, Write(feq_sign))
        self.wait()
        self.play(sec_gp.shift, 0.001 * UP)
        result[2].next_to(feq_sign, direction = RIGHT, buff = 0.25)
        result[1].next_to(result[2], direction = RIGHT, buff = 0.25)
        result[0].next_to(result[1], direction = RIGHT, buff = 0.25)
        self.remove(sec_gp)
        self.play(ReplacementTransform(sec_gp.copy(), result[2]))
        self.wait()
        self.play(Write(result[1]))
        self.wait()
        self.remove(tri_gp)
        self.play(ReplacementTransform(tri_gp.copy(), result[0]))
        self.wait()
        eq_gp = VGroup(integ, feq_sign, result)
        self.play(eq_gp.move_to, ORIGIN, self.camera_frame.move_to, ORIGIN)
        self.wait()
        self.play(FadeOut(feq_sign), FadeOut(result), integ.move_to, ORIGIN)
        self.wait()
        examp = VGroup(
            TexMobject("\\int\\limits_{0}^{1}x^{-x}\\,dx"),
            TexMobject("\\int\\limits_{0}^{b}\\frac{dx}{\\sqrt{a^2-x^2}}"),
            TexMobject("\\int\\limits_{0}^{2\\pi}e^{x\\cos\\theta+y\\sin\\theta}\\,d\\theta"),
            TexMobject("\\int\\limits_{0}^{b}\\frac{dx}{\\sqrt{a^2+x^2}}"),
            TexMobject("\\int\\ln(x^2+a^2)\\,dx"),
            TexMobject("\\int\\limits_{0}^{b}\\sqrt{a^2+x^2}\\,dx")
        )
        for ex in examp:
            self.play(Transform(integ, ex))
        self.wait()
        hyp = self.get_graph(lambda t: np.sqrt(t ** 2 + 16), x_min = -5, x_max = 5, color = BLUE)
        self.play(
            self.camera_frame.move_to, 2 * UP,
            integ.shift, UP,
            Write(hyp)
        )
        self.wait(5)

if __name__ == "__main__":
    module_name = os.path.basename(__file__)
    #command_A = "manim -p  --video_dir ~/Downloads/  "
    clear_cmd = "cls"
    command_A = "python -m manim "
    #command_B = module_name + " " + "LearningScene" + " -pl"
    command_B = module_name + " -p"
    os.system(clear_cmd)
    os.system(command_A + command_B)