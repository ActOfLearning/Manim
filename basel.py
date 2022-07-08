from cmath import sqrt
from manim import *
from numpy import left_shift, right_shift, tri
from scipy import rand

class CreateLight(Scene):
    def construct(self):
        bub = SVGMobject("/Users/muthuveerappanr/Downloads/Manim/BaselProblem/lighthouse-svgrepo-com").set_color(WHITE)
        def light(n = 33, col = YELLOW, fac = 3):
            lights = VMobject()
            for k in range(n):
                temp = Circle(radius = fac * (k + 1) / n, color = col, fill_opacity = (1 - (k / n) ** (1 / 20)), stroke_width = 0.001)
                lights.add(temp)
            lights.add(bub.copy().fade(3 / 8).scale(1 / 2).shift(DOWN / 2 + UP / 8))
            return lights
        temp = light()
        self.play(FadeIn(temp), run_time = 3)
        self.wait()

class Introduction(MovingCameraScene):
    def construct(self):
        bub = SVGMobject("/Users/muthuveerappanr/Downloads/Manim/BaselProblem/lighthouse-svgrepo-com").set_color(WHITE)
        def lighthouse(n = 33, col = YELLOW, fac = 3):
            lights = VMobject()
            for k in range(n):
                temp = Circle(radius = fac * (k + 1) / n, color = col, fill_opacity = (1 - (k / n) ** (1 / 20)), stroke_width = 0.001)
                lights.add(temp)
            lights.add(bub.copy().fade(3 / 8).scale(1 / 2).shift(DOWN / 2 + UP / 8))
            return lights
        temp_lighthosue = lighthouse().scale(1 / 2)
        self.camera.frame.shift(1.5 * UP)

        lake_title = MathTex("\\text{Lake}").scale(2).shift(1.5 * UP)
        intro_circle = Circle(radius = 4, color = BLUE, fill_color = BLUE, fill_opacity = 1 / 8, stroke_width = 4).shift(1.5 * UP)
        self.play(
            Write(lake_title),
            FadeIn(intro_circle),
            run_time = 1.5
        )
        self.wait()
        first_circle = Circle(radius = 1.5, color = BLUE, fill_color = BLUE, fill_opacity = 1 / 8)
        first_lighthouse = temp_lighthosue.copy().fade(0 / 1).move_to(first_circle.get_top())
        observer_point = Dot(radius = 1 / 8, color = BLUE).shift(1.5 * DOWN)
        observer_text = MathTex("\\text{Observer}").next_to(observer_point, DOWN)
        self.play(
            ReplacementTransform(intro_circle, first_circle), 
            FadeOut(lake_title),
            run_time = 1.5
        )
        self.play(
            FadeIn(first_lighthouse),
            Write(observer_point),
            Write(observer_text)
        )
        self.wait()

        def get_future_generations(gen):
            circle_radius = 1.5 * 2 ** gen
            future_circle = Circle(radius = circle_radius, color = BLUE, fill_color = BLUE, fill_opacity = 1 / 8).shift(circle_radius * UP + 1.5 * DOWN)
            circle_center = future_circle.get_center()
            generation_group = VMobject()
            generation_group.add(future_circle)
            temp_dot = Dot().move_to(future_circle.get_bottom())
            temp_dot.rotate(PI / 2 ** gen, about_point = circle_center)
            temp_lighthosue = lighthouse().scale(1 / 1)
            for k in range(2 ** gen):
                future_lighthouse = temp_lighthosue.copy().move_to(temp_dot.get_center())
                temp_dot.rotate(2 * PI / 2 ** gen, about_point = circle_center)
                generation_group.add(future_lighthouse)
            return generation_group
        
        generations = VGroup(VGroup(first_circle, first_lighthouse))
        for k in range(1, 7):
            generations.add(get_future_generations(k))
        self.play(AnimationGroup(FadeOut(generations[0]), FadeIn(generations[1])), run_time = 2)
        for k in range(2, 7):
            generations[k][0].set_style(stroke_width = 8)
        self.play(
            AnimationGroup(FadeOut(generations[1]), FadeIn(generations[2])),
            self.camera.frame.animate.scale(2).shift(3 * UP),
            run_time = 2
        )
        for k in range(2, 7):
            generations[k][0].set_style(stroke_width = 16)
        self.play(
            AnimationGroup(FadeOut(generations[2]), FadeIn(generations[3])),
            self.camera.frame.animate.scale(2).shift(6 * UP),
            run_time = 2
        )
        self.play(AnimationGroup(FadeOut(generations[3]), FadeIn(generations[4])), run_time = 2)
        self.play(AnimationGroup(FadeOut(generations[4]), FadeIn(generations[5])), run_time = 2)
        self.play(AnimationGroup(FadeOut(generations[5]), FadeIn(generations[6])), run_time = 2)
        self.wait()
        basel_text = MathTex("\\frac{1}{1^2}+\\frac{1}{2^2}+\\frac{1}{3^2}+\\cdots=\\frac{\\pi^2}{6}").scale(6).shift(18 * UP)
        basel_text.set_background_stroke(width = 5 * 6)
        wallis_text = MathTex("\\frac{2}{1}\\cdot\\frac{2}{3}\\cdot\\frac{4}{3}\\cdot\\frac{4}{5}\\cdots=\\frac{\\pi}{2}").scale(6)
        wallis_text.set_background_stroke(width = 30)
        wallis_text.next_to(basel_text, DOWN, buff = 2)
        self.play(Write(basel_text))
        self.play(Write(wallis_text))
        self.wait()

        counterparts = VGroup(
            MathTex("\\frac{1}{1^2}+\\frac{1}{4^2}+\\frac{1}{6^2}+\\frac{1}{9^2}+\\cdots"),
            MathTex("\\frac{1}{1}-\\frac{1}{4}+\\frac{1}{6}-\\frac{1}{9}+\\cdots")
        )
        counterparts.arrange(DOWN, buff = 1)
        counterparts.shift(14 * UP)
        counterparts.scale(6)
        counterparts[0].set_color(YELLOW)
        counterparts[1].set_color(TEAL)
        self.play(FadeOut(VGroup(basel_text, wallis_text), shift = UP))
        self.play(Write(counterparts[0]))
        self.play(Write(counterparts[1]))
        self.wait()

class Introduction1(MovingCameraScene):
    def construct(self):
        bub = SVGMobject("/Users/muthuveerappanr/Downloads/Manim/BaselProblem/lighthouse-svgrepo-com").set_color(WHITE)
        first_radius = 1.5
        self.camera.frame.scale(1 + 1 / 8).shift(0.5 * UP)
        def lighthouse(rad = 10.0, col = YELLOW):
            levels = int(10 * rad)
            lights = VMobject()
            dr = rad / levels
            for r in np.arange(0, rad, dr):
                alpha = 1.0 * 1.0 / (r + 1) ** 2
                ann = Annulus(
                    inner_radius = r,
                    outer_radius = r + dr,
                    color = col,
                    fill_opacity = alpha
                )
                lights.add(ann)
            lights.add(bub.copy().fade(3 / 8).scale(3 / 8).shift(DOWN / 2 + UP / 8))
            return lights
        temp_lighthouse = lighthouse().shift(0.5 * UP)
        lake_title = MathTex("\\text{Lake}").scale(2)
        intro_circle = Circle(radius = 4, color = BLUE, fill_color = BLUE, fill_opacity = 1 / 8, stroke_width = 4)
        self.play(
            Write(lake_title),
            FadeIn(intro_circle),
            run_time = 1.5
        )
        self.wait()

        observer_point = Dot(color = WHITE).shift(3.5 * DOWN)
        observer_title = MathTex("\\text{Observer}").next_to(observer_point, DOWN, buff = 1 / 16)
        self.play(
            FadeOut(lake_title),
            intro_circle.animate.scale(1 / 2).shift(1.5 * DOWN)
        )
        self.play(
            FadeIn(temp_lighthouse),
            Write(observer_point),
            Write(observer_title)
        )
        self.wait()

        def get_future_generations(gen):
            circle_radius = 2 * 2 ** gen
            future_circle = Circle(radius = circle_radius, color = BLUE, fill_color = BLUE, fill_opacity = 1 / 8).shift(circle_radius * UP + 3.5 * DOWN)
            circle_center = future_circle.get_center()
            generation_group = VMobject()
            generation_group.add(future_circle)
            temp_dot = Dot().move_to(future_circle.get_bottom())
            temp_dot.rotate(PI / 2 ** gen, about_point = circle_center)
            if gen > 3:
                temp_dot.rotate(-2 * PI / 2 ** gen - 5 * 2 * PI / 2 ** gen, about_point = circle_center)
            temp_lighthosue = lighthouse()
            for k in range(min(10, 2 ** gen)):
                future_lighthouse = temp_lighthosue.copy().move_to(temp_dot.get_center())
                temp_dot.rotate(2 * PI / 2 ** gen, about_point = circle_center)
                generation_group.add(future_lighthouse)
            return generation_group
        
        generations = VGroup(VGroup(intro_circle, temp_lighthouse))
        generations.add(get_future_generations(1))
        self.play(Write(generations[1][0]))
        self.play(
            FadeOut(generations[0]),
            GrowFromCenter(generations[1][1:])
        )
        self.wait()

        self.play(generations[1].animate.scale(1 / 2, about_point = 3.5 * DOWN))

        generations.add(get_future_generations(2).scale(1 / 2, about_point = 3.5 * DOWN))
        self.play(Write(generations[2][0]))
        self.play(
            FadeOut(generations[1]),
            GrowFromCenter(generations[2][1:])
        )
        self.wait()

        self.play(generations[2].animate.scale(1 / 2, about_point = 3.5 * DOWN))
        self.wait()

        generations.add(get_future_generations(3).scale(1 / 4, about_point = 3.5 * DOWN))
        self.play(
            FadeOut(generations[2]),
            FadeIn(generations[3])
        )
        self.wait()

        for k in range(5):
            generations.add(get_future_generations(k + 4).scale(1 / 4, about_point = 3.5 * DOWN))
            self.play(
                FadeOut(generations[3 + k]),
                FadeIn(generations[k + 4])
            )
            self.wait()
        self.wait()

class Introduction2(MovingCameraScene):
    def construct(self):
        bub = SVGMobject("/Users/muthuveerappanr/Downloads/Manim/BaselProblem/lighthouse-svgrepo-com").set_color(WHITE)
        def lighthouse(rad = 10.0, col = YELLOW):
            levels = int(10 * rad)
            lights = VMobject()
            dr = rad / levels
            for r in np.arange(0, rad, dr):
                alpha = 1.0 * 1.0 / (r + 1) ** 2
                ann = Annulus(
                    inner_radius = r,
                    outer_radius = r + dr,
                    color = col,
                    fill_opacity = alpha
                )
                lights.add(ann)
            lights.add(bub.copy().fade(3 / 8).scale(3 / 8).shift(DOWN / 2 + UP / 8))
            return lights
        temp_lighthouse = lighthouse().shift(6 * RIGHT)
        temp_dot = Dot().shift(6 * RIGHT)
        lighthouses = VGroup()
        n = 10
        for k in range(n):
            lighthouses.add(temp_lighthouse.copy().move_to(temp_dot.get_center()))
            temp_dot.rotate_about_origin(2 * PI / n)
        self.play(LaggedStart(*[FadeIn(obj) for obj in lighthouses]), run_time = 5)
        self.wait()

        basel_problem = MathTex("\\frac{1}{1^2}+\\frac{1}{2^2}+\\frac{1}{3^2}+\\cdots", "=", "\\frac{\\pi^2}{6}")
        wallis_product = MathTex("\\frac{2}{1}\\cdot\\frac{2}{3}\\cdot\\frac{4}{3}\\cdot\\frac{4}{5}\\cdots", "=", "\\frac{\\pi}{2}")
        basel_problem.set_background_stroke(width = 5)
        wallis_product.set_background_stroke(width = 5)

        self.play(Write(basel_problem))
        self.wait()

        other_lighthouse = lighthouse(col = BLUE)
        other_lighthouses = VGroup()
        other_dot = Dot().shift(6 * RIGHT).rotate_about_origin(PI / n)
        for k in range(n):
            other_lighthouses.add(other_lighthouse.copy().move_to(other_dot.get_center()))
            other_dot.rotate_about_origin(2 * PI / n)
        anims = [ReplacementTransform(obj1, obj2) for obj1, obj2 in zip(lighthouses, other_lighthouses)]
        wallis_product.shift(1.5 * DOWN)
        self.play(
            AnimationGroup(*anims, run_time = 3, lag_ratio = 0.1),
            basel_problem.animate.shift(1.5 * UP),
            FadeIn(wallis_product),
        )
        self.wait()

        basel_like1 = MathTex("\\frac{1}{1^2}+\\frac{1}{4^2}+\\frac{1}{6^2}+\\frac{1}{9^2}+\\cdots").shift(1.5 * UP)
        linear_basel = MathTex("\\frac{1}{1}-\\frac{1}{4}+\\frac{1}{6}-\\frac{1}{9}+\\cdots").shift(1.5 * DOWN)
        self.play(
            AnimationGroup(FadeOut(basel_problem, shift = RIGHT), Write(basel_like1), lag_ratio = 0.25),
            AnimationGroup(FadeOut(wallis_product, shift = RIGHT), Write(linear_basel), lag_ratio = 0.25)
        )
        self.wait()

class OtherPointOnCircle(MovingCameraScene):
    def construct(self):
        bub = SVGMobject("/Users/muthuveerappanr/Downloads/Manim/BaselProblem/lighthouse-svgrepo-com").set_color(WHITE)
        def lighthouse(rad = 10.0, col = YELLOW):
            levels = int(10 * rad)
            lights = VMobject()
            dr = rad / levels
            for r in np.arange(0, rad, dr):
                alpha = 1.0 * 1.0 / (r + 1) ** 2
                ann = Annulus(
                    inner_radius = r,
                    outer_radius = r + dr,
                    color = col,
                    fill_opacity = alpha
                )
                lights.add(ann)
            lights.add(bub.copy().fade(3 / 8).scale(3 / 8).shift(DOWN / 2 + UP / 8))
            return lights
        single_lighthouse = lighthouse()
        lake = Circle(radius = 3, color = BLUE, fill_color = BLUE, fill_opacity = 0.25)
        lighthouse_dot = Dot(fill_opacity = 0).shift(lake.get_top())
        single_lighthouse.move_to(lighthouse_dot.get_center())
        observer = Dot().move_to(lake.get_bottom())
        self.play(FadeIn(lake))
        self.play(FadeIn(observer), FadeIn(single_lighthouse))
        single_lighthouse.add_updater(lambda obj: obj.move_to(lighthouse_dot))
        self.wait()

        self.play(
            Rotate(lighthouse_dot, about_point = ORIGIN, angle = PI / 4),
            run_time = 2
        )
        self.play(
            Rotate(lighthouse_dot, about_point = ORIGIN, angle = -PI / 2),
            run_time = 2
        )
        self.play(
            Rotate(lighthouse_dot, about_point = ORIGIN, angle = 2 * PI / 5 - PI / 4 - PI / 2),
            run_time = 2
        )

        diameter_line = Line(lake.get_top(), lake.get_bottom())
        connecting_line = Line(lake.get_top(), lighthouse_dot.get_center())
        ang = Angle(diameter_line, connecting_line, radius = 0.5)
        pibyfive = MathTex("\\frac{\\pi}{5}").next_to(ang, direction = DOWN).shift(RIGHT / 8)
        pibyfive.set_color(YELLOW)
        pibyfive.set_background_stroke(width = 5)
        self.play(Write(diameter_line), Write(connecting_line))
        self.play(Write(ang), Write(pibyfive))
        self.wait()

        brace = Brace(diameter_line, LEFT)
        fivebypi = MathTex("\\frac{5}{\\pi}").next_to(brace, LEFT).set_background_stroke(width = 5)
        self.play(Write(brace), Write(fivebypi))
        self.wait()

        brightness_circle = Circle(radius = 0.75, color = WHITE, fill_color = YELLOW, fill_opacity = 0.5)
        brighness_measure = MathTex("\\frac{1}{r^2}").set_background_stroke(width = 5)
        brightness_group = VGroup(brightness_circle, brighness_measure)
        brightness_group.move_to(lake.get_bottom() + 1 * DOWN / 4).shift(2.5 * LEFT)
        self.play(
            Write(brightness_group),
        )
        self.wait()

        observer_lighthouse_line = Line(observer, lighthouse_dot)
        right_angle = Square(color = WHITE, fill_opacity = 0.75).scale(1 / 4)
        right_angle.move_to(lighthouse_dot)
        right_angle.shift(right_angle.get_center() - right_angle.get_critical_point(DL))
        right_angle.rotate(PI / 5 + PI / 2, about_point = right_angle.get_critical_point(DL))
        self.play(
            Write(observer_lighthouse_line),
            Write(right_angle)
        )
        self.wait()

        rtext = MathTex("r")
        rtext.move_to(observer_lighthouse_line.get_center() + 3 * UP / 8)
        rtext.set_background_stroke(width = 5)
        requals = MathTex("r=\\frac{5}{\\pi}\\sin\\left(\\frac{\\pi}{5}\\right)")
        requals.to_corner(DR)
        self.play(
            Write(rtext),
            Write(requals)
        )
        self.wait()

        self.play(
            FadeOut(right_angle),
            FadeOut(rtext),
            FadeOut(requals),
            FadeOut(observer_lighthouse_line)
        )
        self.wait()

        center_dot = Dot(lake.get_center())
        radial_line = Line(center_dot, lighthouse_dot)
        center_ang = Angle(diameter_line, radial_line)
        twopibyfive = MathTex("\\frac{2\\pi}{5}").set_color(YELLOW)
        twopibyfive.set_background_stroke(width = 5)
        twopibyfive.next_to(center_ang, DOWN).shift(RIGHT / 4)
        insc = MathTex("\\text{Inscribed angle}", "\\text{theorem}").to_corner(UL)
        insc[-1].shift(2.5 * LEFT + DOWN / 2)
        self.play(
            Write(center_dot),
            Write(radial_line),
            Write(center_ang),
            Write(twopibyfive),
            Write(insc)
        )
        self.wait()

        temp_arc = Arc(radius = 3, start_angle = 0, angle = 2 * PI / 5, color = WHITE)
        arcbrace = ArcBrace(temp_arc).rotate(-PI / 2, about_point = lake.get_center())
        arclength = MathTex("\\frac{1}{2}\\cdot\\frac{5}{\\pi}\\cdot \\frac{2\\pi}{5}=1")
        arclength.next_to(arcbrace)
        self.play(
            Write(arcbrace),
            Write(arclength)
        )
        self.wait()

class RecapALittle(MovingCameraScene):
    def construct(self):
        bub = SVGMobject("/Users/muthuveerappanr/Downloads/Manim/BaselProblem/lighthouse-svgrepo-com").set_color(WHITE)
        def lighthouse(rad = 10.0, col = YELLOW):
            levels = int(10 * rad)
            lights = VMobject()
            dr = rad / levels
            for r in np.arange(0, rad, dr):
                alpha = 1.0 * 1.0 / (r + 1) ** 2
                ann = Annulus(
                    inner_radius = r,
                    outer_radius = r + dr,
                    color = col,
                    fill_opacity = alpha
                )
                lights.add(ann)
            lights.add(bub.copy().fade(3 / 8).scale(3 / 8).shift(DOWN / 2 + UP / 8))
            return lights
        co_axes = VGroup(Line(20 * LEFT, 20 * RIGHT), Line(10 * DOWN, 10 * UP))
        co_axes.set_color(GREY)
        self.play(Write(co_axes))
        self.play(self.camera.frame.animate.move_to(2.75 * UP + 5.5 * RIGHT))
        self.wait()

        x_val, y_val = 8, 5
        px_val, py_val = 2.247, 3.595
        x_pos, y_pos = x_val * RIGHT, y_val * UP
        x_lighthouse = lighthouse()
        y_lighthouse = lighthouse()
        single_lighthouse = lighthouse()
        xy_line = Line(x_pos, y_pos)
        x_line, y_line = Line(ORIGIN, x_pos).set_color(GREEN), Line(ORIGIN, y_pos).set_color(RED)
        p_line = Line(ORIGIN, px_val * RIGHT + py_val * UP, color = GREY)
        little_a, little_b = MathTex("a").set_color(GREEN), MathTex("b").set_color(RED)
        little_h = MathTex("h").set_color(YELLOW)
        little_a.next_to(x_line, DOWN, buff = 1 / 4)
        little_b.next_to(y_line, LEFT, buff = 1 / 4)
        little_h.move_to(p_line.get_center()).shift(RIGHT / 2)
        for obj in [single_lighthouse, x_lighthouse, y_lighthouse]:
            obj.shift(px_val * RIGHT + py_val * UP)
        self.play(FadeIn(single_lighthouse), Write(p_line))
        self.wait()
        self.play(Write(xy_line))
        self.play(
            Write(x_line),
            Write(y_line),
            Write(little_a),
            Write(little_b),
            Write(little_h)
        )
        self.wait()
        self.add(x_lighthouse, y_lighthouse)
        self.play(
            #FadeOut(single_lighthouse),
            x_lighthouse.animate.move_to(x_pos * RIGHT),
            y_lighthouse.animate.move_to(y_pos * UP),
        )
        self.wait()

        inv_pyth = MathTex("\\frac{1}{h^2}", "=", "\\frac{1}{a^2}", "+", "\\frac{1}{b^2}")
        title = MathTex("\\text{Inverse Pythagorean Theorem}").to_edge(UP).shift(2.75 * UP + 5.5 * RIGHT)
        title_line = Underline(title)
        title_group = VGroup(title, title_line)
        inv_pyth.scale(1.25)
        inv_pyth[0][-2].set_color(YELLOW)
        inv_pyth[2][-2].set_color(GREEN)
        inv_pyth[-1][-2].set_color(RED)
        inv_pyth.shift(8 * RIGHT + 3.5 * UP)
        self.play(Write(inv_pyth), run_time = 2)
        self.play(Write(title_group))
        self.wait()

class LikeInOtherVideo(MovingCameraScene):
    def construct(self):
        bub = SVGMobject("/Users/muthuveerappanr/Downloads/Manim/BaselProblem/lighthouse-svgrepo-com").set_color(WHITE)
        def lighthouse(rad = 10.0, col = YELLOW):
            levels = int(10 * rad)
            lights = VMobject()
            dr = rad / levels
            for r in np.arange(0, rad, dr):
                alpha = 1.0 * 1.0 / (r + 1) ** 2
                ann = Annulus(
                    inner_radius = r,
                    outer_radius = r + dr,
                    color = col,
                    fill_opacity = alpha
                )
                lights.add(ann)
            lights.add(bub.copy().fade(3 / 8).scale(3 / 8).shift(DOWN / 2 + UP / 8))
            return lights
        single_lighthouse = lighthouse()
        lighthouse_dot = Dot(fill_opacity = 0)
        lake = Circle(radius = 3.5, color = BLUE, fill_color = BLUE, fill_opacity = 1 / 8)
        lake1 = lake.copy()
        lake2 = lake.copy()
        observer = Dot().move_to(lake.get_bottom())
        lighthouse_dot.move_to(lake.get_bottom())
        lighthouse_dot.rotate_about_origin(2 * PI / 5)
        single_lighthouse.move_to(lighthouse_dot.get_center())
        self.play(Write(lake))
        self.play(FadeIn(single_lighthouse), Write(observer))
        self.wait()
        self.play(VGroup(lake, single_lighthouse).animate.scale(1 / 2, about_point = observer.get_center()))
        self.wait()
        self.play(Write(lake1))
        centerdot1 = Dot().move_to(lake1.get_center())
        diameter1 = Line(lake1.get_top(), lake1.get_bottom())
        diameter1.rotate_about_origin(PI / 5)
        self.play(Write(diameter1), Write(centerdot1))
        self.wait()
        conn_line1, conn_line2 = Line(observer.get_center(), diameter1.get_start()), Line(observer.get_center(), diameter1.get_end())
        right_angle = Square(color = WHITE, fill_opacity = 0.75).scale(1 / 4)
        right_angle.move_to(observer.get_center() + right_angle.get_center() - right_angle.get_critical_point(DL))
        right_angle.rotate(PI / 10, about_point = right_angle.get_critical_point(DL))
        self.play(
            Write(conn_line1),
            Write(conn_line2)
        )
        self.play(FadeIn(right_angle))
        lighthouse1, lighthouse2 = single_lighthouse.copy(), single_lighthouse.copy()
        self.play(
            lighthouse1.animate.move_to(diameter1.get_start()),
            lighthouse2.animate.move_to(diameter1.get_end()),
        )
        self.wait()
        vertical_radius1 = Line(lake1.get_center(), lake1.get_bottom(), color = BLUE)
        vertical_radius1_brace = Brace(vertical_radius1, LEFT)
        diameter = Line(lake.get_top(), lake.get_bottom(), color = BLUE)
        fivebypi = MathTex("\\frac{5}{\\pi}").set_background_stroke(width = 5)
        fivebypi.next_to(vertical_radius1_brace, LEFT)
        pibyfive_angle = Angle(diameter, Line(lake1.get_center(), diameter1.get_end(), color = RED))
        pibyfive = MathTex("\\frac{\\pi}{5}").next_to(lake1.get_center())
        pibyfive.set_background_stroke(width = 5)
        half_conn_line1 = Line(lake1.get_center(), diameter1.get_end())
        self.add(half_conn_line1)
        self.play(
            Write(vertical_radius1),
            FadeOut(right_angle),
            FadeOut(conn_line1),
            FadeOut(conn_line2),
            FadeOut(diameter1)
        )
        self.wait()
        self.play(
            Write(vertical_radius1_brace),
            Write(fivebypi),
        )
        self.wait()
        self.play(
            Write(pibyfive_angle),
            Write(pibyfive)
        )
        pibyten_arc = Arc(radius = 3.5, start_angle = 0, angle = PI / 5)
        pibyten_arc_brace = ArcBrace(pibyten_arc).rotate(-PI / 2, about_point = lake1.get_center())
        one_unit = MathTex("1\\text{ unit}").move_to(pibyten_arc_brace.get_center())
        one_unit.shift(DOWN / 4 + 1.5 * RIGHT)
        one_unit.set_background_stroke(width = 5)
        self.play(
            Write(pibyten_arc_brace),
            Write(one_unit)
        )
        self.wait()
        circumf = MathTex("\\text{Circumference}=10\\text{ units}").set_background_stroke(width = 10)
        circumf.move_to(lake1.get_center()).shift(UP)
        self.play(
            Write(diameter1),
            FadeOut(lake),
            FadeOut(single_lighthouse),
            FadeOut(pibyfive),
            FadeOut(pibyfive_angle),
            FadeOut(pibyten_arc),
            FadeOut(pibyten_arc_brace),
            FadeOut(one_unit),
        )
        self.play(Write(circumf), run_time = 2)
        self.wait()
        centerdot2 = Dot(lake2.get_center())
        self.play(
            FadeOut(vertical_radius1),
            FadeOut(vertical_radius1_brace),
            FadeOut(fivebypi),
            FadeOut(half_conn_line1),
            FadeOut(circumf),
            FadeOut(diameter1),
            FadeOut(centerdot1),
            Write(lake2),
            VGroup(lake1, lighthouse1, lighthouse2).animate.scale(1 / 2, about_point = lake.get_bottom())
        )
        self.play(Write(centerdot2))
        self.wait()

class RepeatingTheProcess(MovingCameraScene):
    def construct(self):
        bub = SVGMobject("/Users/muthuveerappanr/Downloads/Manim/BaselProblem/lighthouse-svgrepo-com").set_color(WHITE)
        def lighthouse(rad = 10.0, col = YELLOW, fac = 1.0, lfac = 10):
            levels = int(lfac * rad)
            lights = VMobject()
            dr = rad / levels
            for r in np.arange(0, rad, dr):
                alpha = fac * 1.0 / (r + 1) ** 2
                ann = Annulus(
                    inner_radius = r,
                    outer_radius = r + dr,
                    color = col,
                    fill_opacity = alpha
                )
                lights.add(ann)
            lights.add(bub.copy().fade(3 / 8).scale(3 / 8).shift(DOWN / 2 + UP / 8))
            return lights
        lake = Circle(radius = 3.5, color = BLUE, fill_color = BLUE, fill_opacity = 1 / 8)
        observer = Dot().move_to(lake.get_bottom())
        def get_future_generations(gen):
            base_circle = Circle(radius = 3.5, color = BLUE, fill_color = BLUE, fill_opacity = 1 / 8)
            circle_center = base_circle.get_center()
            generation_group = VMobject()
            generation_group.add(base_circle)
            temp_dot = Dot().move_to(base_circle.get_bottom())
            temp_dot.rotate(2 * PI / 5 / 2 ** gen, about_point = circle_center)
            temp_lighthosue = lighthouse().scale(1 / 2)
            for k in range(2 ** gen):
                future_lighthouse = temp_lighthosue.copy().move_to(temp_dot.get_center())
                temp_dot.rotate(2 * PI / 2 ** gen, about_point = circle_center)
                generation_group.add(future_lighthouse)
            return generation_group
        def get_lines(gen):
            res = VGroup()
            line = Line(3.5 * DOWN, 3.5 * UP, stroke_width = 1 / 1)
            line.rotate_about_origin(2 * PI / 5 / 2 ** gen)
            for k in range(2 ** gen):
                res.add(line.copy().rotate_about_origin(k * 2 * PI / 2 ** gen))
            return res
        first_gen = get_future_generations(1)
        self.play(FadeIn(first_gen), Write(observer))
        self.play(first_gen.animate.scale(1 / 2, about_point = first_gen[0].get_bottom()))
        self.wait()

        second_gen = get_future_generations(2)
        second_gen_center = Dot().move_to(second_gen[0].get_center())
        self.play(Write(second_gen[0]), Write(second_gen_center))
        self.wait()
        lines = VGroup()
        for k in range(2):
            new_lighthouses = VGroup(second_gen[k + 1], second_gen[k + 1 + 2])
            new_lighthouses.save_state()
            lines.become(get_lines(2))
            self.play(Write(lines[k]))
            for obj in new_lighthouses:
                obj.move_to(first_gen[k + 1])
            self.add(new_lighthouses)
            self.play(Restore(new_lighthouses))
            self.wait()
        second_circum = MathTex("\\text{Circumference: }20\\text{ units}").set_background_stroke(width = 5)
        self.play(Write(second_circum), run_time = 2)
        self.wait()
        self.play(
            FadeOut(lines),
            FadeOut(first_gen),
            FadeOut(second_gen_center),
            FadeOut(second_circum, shift = DOWN),
            second_gen.animate.scale(1 / 2, about_point = second_gen[0].get_bottom())
        )
        self.wait()

        third_gen = get_future_generations(3)
        third_gen_center = Dot().move_to(third_gen[0].get_center())
        self.play(Write(third_gen[0]), Write(third_gen_center))
        self.wait()
        for k in range(4):
            new_lighthouses = VGroup(third_gen[k + 1], third_gen[k + 1 + 4])
            new_lighthouses.save_state()
            lines.become(get_lines(3))
            self.play(Write(lines[k]))
            for obj in new_lighthouses:
                obj.move_to(second_gen[k + 1])
            self.add(new_lighthouses)
            self.play(Restore(new_lighthouses))
            self.wait()
        third_circum = MathTex("\\text{Circumference: }40\\text{ units}").set_background_stroke(width = 5)
        self.play(Write(third_circum))
        self.wait()
        self.play(
            FadeOut(lines),
            FadeOut(second_gen),
            FadeOut(third_gen_center),
            FadeOut(third_circum, shift = DOWN),
            #third_gen.animate.scale(1 / 2, about_point = third_gen[0].get_bottom())
        )
        self.wait()
        def get_partial_generations(gen):
            base_circle = Circle(radius = 3.5 * 2 ** (gen - 3), color = BLUE, fill_color = BLUE, fill_opacity = 1 / 8)
            base_circle.move_to(3.5 * DOWN + base_circle.get_center() - base_circle.get_bottom())
            circle_center = base_circle.get_center()
            generation_group = VMobject()
            generation_group.add(base_circle)
            temp_dot = Dot().move_to(base_circle.get_bottom())
            temp_dot.rotate(2 * PI / 5 / 2 ** gen, about_point = circle_center)
            temp_dot.rotate(-4 * 2 * PI / 2 ** gen, about_point = circle_center)
            temp_lighthosue = lighthouse().scale(1 / 2)
            for k in range(9):
                future_lighthouse = temp_lighthosue.copy().move_to(temp_dot.get_center())
                temp_dot.rotate(2 * PI / 2 ** gen, about_point = circle_center)
                generation_group.add(future_lighthouse)
            return generation_group
        fourth_partial_gen = get_partial_generations(4)
        self.play(
            FadeOut(third_gen),
            FadeIn(fourth_partial_gen)
        )
        self.wait()
        fifth_partial_gen = get_partial_generations(5)
        self.play(
            FadeOut(fourth_partial_gen),
            FadeIn(fifth_partial_gen)
        )
        self.wait()
        sixth_partial_gen = get_partial_generations(6)
        self.play(
            FadeOut(fifth_partial_gen),
            FadeIn(sixth_partial_gen)
        )
        self.wait()
        seventh_partial_gen = get_partial_generations(7)
        self.play(
            FadeOut(sixth_partial_gen),
            FadeIn(seventh_partial_gen)
        )
        self.wait()
        eigth_partial_gen = get_partial_generations(8)
        self.play(
            FadeOut(seventh_partial_gen),
            FadeIn(eigth_partial_gen)
        )
        self.wait()
        circle_line = Line(10 * LEFT, 10 * RIGHT, color = BLUE)
        line_lighthouses = VGroup(circle_line)
        for k in range(-9, 16, 5):
            temp_lighthouse = lighthouse().scale(1 / 2)
            #temp_lighthouse.move_to(observer.get_center())
            temp_lighthouse.shift(3.5 * PI / 20 * k * RIGHT)
            line_lighthouses.add(temp_lighthouse)
        line_lighthouses.shift(3.5 * DOWN)
        ticks = VGroup()
        tick = Line(UP / 8, DOWN / 8, color = GREEN)
        for k in range(-15, 23):
            j = k % 5
            if j == 1:
                continue
            ticks.add(tick.copy().shift(3.5 * DOWN + 3.5 * PI / 20 * k * RIGHT))
        self.play(
            FadeOut(eigth_partial_gen),
            FadeIn(line_lighthouses)
        )
        self.wait()
        onecommazero = MathTex("(1,0)").set_background_stroke(width = 5)
        onecommazero.shift(3.5 * DOWN + 3.5 * PI / 20 * 1 * RIGHT + UR / 4)
        self.play(Write(onecommazero))
        self.wait()
        self.play(ShowIncreasingSubsets(ticks), run_time = 3)
        self.wait()
        zero_gen = get_future_generations(0)
        self.play(FadeIn(zero_gen))
        self.wait()
        triangle_group = VGroup(
            Line(zero_gen[0].get_top(), zero_gen[0].get_bottom()),
            Line(zero_gen[0].get_bottom(), zero_gen[1][0].get_center()),
            Line(zero_gen[1][0].get_center(), zero_gen[0].get_top()),
        )
        right_angle = Square(color = WHITE, fill_opacity = 0.75).scale(1 / 4)
        right_angle.move_to(zero_gen[1][0].get_center() + right_angle.get_center() - right_angle.get_critical_point(DR))
        right_angle.rotate(PI / 5, about_point = right_angle.get_critical_point(DR))
        self.play(Write(triangle_group), run_time = 2)
        self.wait()
        brace = Brace(triangle_group[0], LEFT)
        fivebypi = MathTex("\\frac{5}{\\pi}").set_background_stroke(width = 5)
        fivebypi.next_to(brace, LEFT)
        pibyfive = MathTex("\\frac{\\pi}{5}").set_background_stroke(width = 5)
        angle = Angle(triangle_group[0], Line(triangle_group[-1].get_end(), triangle_group[-1].get_start()))
        pibyfive.next_to(angle, DOWN).shift(RIGHT / 8)
        self.play(
            FadeIn(brace),
            FadeIn(fivebypi),
            FadeIn(pibyfive),
            FadeIn(angle),
            FadeIn(right_angle)
        )
        self.wait()
        sine_length = MathTex("\\frac{5}{\\pi}\\sin\\left(\\frac{\\pi}{5}\\right)").set_color(ORANGE).set_background_stroke(width = 5)
        sine_length.rotate_about_origin(PI / 5)
        sine_length.move_to(triangle_group[1].get_center())
        sine_length.shift(LEFT / 2 + UP / 2)
        self.play(
            Write(sine_length),
            FadeOut(onecommazero),
            FadeOut(fivebypi),
            FadeOut(brace),
            FadeOut(pibyfive),
            FadeOut(angle),
            triangle_group[1].animate.set_color(ORANGE)
        )
        self.wait()
        self.play(
            FadeOut(right_angle),
            FadeOut(VGroup(triangle_group[0], triangle_group[-1]))
        )
        self.wait()
        final_eqn = MathTex("\\frac{1}{\\left(\\frac{5}{\\pi}\\right)^2\\sin^2\\left(\\frac{\\pi}{5}\\right)}", "=", "\\cdots+\\frac{1}{(-9)^2}+\\frac{1}{(-4)^2}+\\frac{1}{1^2}+\\frac{1}{6^2}+\\frac{1}{11^2}+\\cdots")
        final_eqn.set_background_stroke(width = 5)
        final_eqn.shift(UP)
        self.play(Write(final_eqn), run_time = 3)
        self.wait()
        final_eqn2 = MathTex("\\frac{(\\pi/5)^2}{\\sin^2\\left(\\frac{\\pi}{5}\\right)}", "=", "\\frac{1}{1^2}+\\frac{1}{4^2}+\\frac{1}{6^2}+\\frac{1}{9^2}+\\frac{1}{11^2}+\\cdots")
        final_eqn2.set_background_stroke(width = 5)
        final_eqn2.shift(UP)
        self.play(
            final_eqn.animate.shift(2 * UP),
            FadeIn(final_eqn2, shift = UP)
        )
        self.wait()
        self.play(FadeOut(VGroup(triangle_group[1], sine_length)))
        self.wait()
        #final_eqn3 = MathTex("\\frac{(\\pi/5)^2}{\\sin^2\\left(\\frac{\\pi}{5}\\right)}", "=", "\\cdots+\\frac{1}{(\\frac{1}{5}-2)^2}+\\frac{1}{(\\frac{1}{5}-1)^2}+\\frac{1}{(\\frac{1}{5}-0)^2}+\\frac{1}{(\\frac{1}{5}+1)^2}+\\frac{1}{(\\frac{1}{5}+2)^2}+\\cdots")
        final_eqn3 = MathTex("\\frac{\\pi^2}{\\sin^2\\left(\\frac{\\pi}{5}\\right)}", "=", "\\sum_{k \\in \\mathbb{Z}}\\frac{1}{(\\frac{1}{5}+k)^2}")
        final_eqn3.set_background_stroke(width = 5)
        final_eqn3.shift(UP)
        self.play(
            final_eqn.animate.shift(2 * UP),
            final_eqn2.animate.shift(2 * UP),
            FadeIn(final_eqn3, shift = UP)
        )
        self.wait()
        #general_eqn = MathTex("\\frac{\\pi^2}{\\sin^2(\\pi x)}", "=", "\\cdots+\\frac{1}{(x-2)^2}+\\frac{1}{(x-1)^2}+\\frac{1}{(x-0)^2}+\\frac{1}{(x+1)^2}+\\frac{1}{(x+2)^2}+\\cdots")
        general_eqn = MathTex("\\frac{\\pi^2}{\\sin^2(\\pi x)}", "=", "\\sum_{k \\in \\mathbb{Z}}\\frac{1}{(x+k)^2}")
        general_eqn.set_background_stroke(width = 5)
        self.play(
            final_eqn.animate.shift(2 * UP),
            final_eqn2.animate.shift(2 * UP),
            final_eqn3.animate.shift(4 * UP),
            Write(general_eqn, run_time = 3)
        )
        self.wait()
        self.play(
            FadeOut(general_eqn),
            FadeOut(ticks),
            FadeOut(zero_gen)
        )
        self.wait()
        new_light = lighthouse(rad = 10, col = BLUE, fac = 2.5, lfac = 50)
        anims = []
        for obj in new_light:
            anims += [FadeIn(obj)]
        self.play(
            AnimationGroup(*anims),
            run_time = 5
        )
        self.wait()

class NewLighthouses(MovingCameraScene):
    def construct(self):
        bub = SVGMobject("/Users/muthuveerappanr/Downloads/Manim/BaselProblem/lighthouse-svgrepo-com").set_color(WHITE)
        def lighthouse(rad = 10.0, col = YELLOW, fac = 1.0, lfac = 10):
            levels = int(lfac * rad)
            lights = VMobject()
            dr = rad / levels
            for r in np.arange(0, rad, dr):
                alpha = fac * 1.0 / (r + 1) ** 2
                ann = Annulus(
                    inner_radius = r,
                    outer_radius = r + dr,
                    color = col,
                    fill_opacity = alpha
                )
                lights.add(ann)
            lights.add(bub.copy().fade(3 / 8).scale(3 / 8).shift(DOWN / 2 + UP / 8))
            return lights
        horizontal = Line(10 * LEFT, 10 * RIGHT)
        observer = Dot()
        self.play(
            Write(horizontal),
            FadeIn(observer)
        )
        self.wait()
        right_light = lighthouse(col = BLUE)
        right_light.shift(5 * RIGHT)
        left_light = lighthouse(col = YELLOW)
        left_light.shift(5 * LEFT)
        dist_r = MathTex("r").shift(2.5 * RIGHT + DOWN / 4)
        self.play(FadeIn(right_light))
        self.play(Write(dist_r))
        self.wait()
        blue_bright_expression = MathTex("\\text{Brightness}", "\\propto", "\\frac{1}{r}").shift(2 * UP)
        blue_bright_expression[0].set_color(BLUE)
        yellow_bright_expression = MathTex("\\text{Brightness}", "\\propto", "\\frac{1}{r^2}").shift(2 * DOWN)
        yellow_bright_expression[0].set_color(YELLOW)
        self.play(Write(blue_bright_expression), run_time = 2)
        self.wait()
        self.play(FadeIn(left_light))
        self.play(Write(yellow_bright_expression), run_time = 2)
        self.wait()
        r_surrounding = SurroundingRectangle(blue_bright_expression[-1])
        rsq_surrounding = SurroundingRectangle(yellow_bright_expression[-1])
        self.play(Create(rsq_surrounding))
        self.play(FadeTransform(rsq_surrounding, r_surrounding))
        self.play(Uncreate(r_surrounding))
        self.wait()
        right_blue_light = lighthouse(col = BLUE)
        right_blue_light.shift(5 * LEFT)
        neg_bright = MathTex("\\text{Brightness}", "\\propto", "\\frac{1}{-r}").shift(5 * LEFT + 2 * DOWN)
        neg_bright[0].set_color(BLUE)
        self.play(
            FadeOut(dist_r),
            blue_bright_expression.animate.shift(5 * RIGHT),
            Transform(yellow_bright_expression, neg_bright),
            FadeOut(left_light),
            FadeIn(right_blue_light)
        )
        self.wait()
        right_lines, left_lines = VGroup(), VGroup()
        for k in range(-20, 20 + 1):
            right_lines.add(Line(5 * RIGHT, k * UP / 40, stroke_width = 1))
            left_lines.add(Line(5 * LEFT, k * UP / 40, stroke_width = 1))
        right_lines.set_color(BLUE)
        left_lines.set_color(BLUE)
        self.play(Create(right_lines), Create(left_lines), run_time = 5)
        self.play(Uncreate(right_lines), Uncreate(left_lines), run_time = 5)
        self.wait()
        self.play(
            right_light.animate.move_to(ORIGIN),
            right_blue_light.animate.move_to(ORIGIN),
            FadeOut(VGroup(yellow_bright_expression, blue_bright_expression))
        )
        self.wait()

class InverseSumTheorem(MovingCameraScene):
    def construct(self):
        bub = SVGMobject("/Users/muthuveerappanr/Downloads/Manim/BaselProblem/lighthouse-svgrepo-com").set_color(WHITE)
        def lighthouse(rad = 10.0, col = YELLOW, fac = 1.0, lfac = 10):
            levels = int(lfac * rad)
            lights = VMobject()
            dr = rad / levels
            for r in np.arange(0, rad, dr):
                alpha = fac * 1.0 / (r + 1) ** 2
                ann = Annulus(
                    inner_radius = r,
                    outer_radius = r + dr,
                    color = col,
                    fill_opacity = alpha
                )
                lights.add(ann)
            lights.add(bub.copy().fade(3 / 8).scale(3 / 8).shift(DOWN / 2 + UP / 8))
            return lights
        horizontal = Line(20 * LEFT, 20 * RIGHT, color = GREY)
        vertical = Line(20 * UP, 20 * DOWN, color = GREY)
        self.play(
            Write(horizontal),
            Write(vertical),
            self.camera.frame.animate.shift(5 * RIGHT + 2.5 * UP)
        )
        self.wait()
        light1dot, light2dot = Dot(fill_opacity = 0), Dot(fill_opacity = 0)
        light1position, light2position = ValueTracker(5), ValueTracker(8)
        light1, light2 = lighthouse(col = BLUE), lighthouse(col = BLUE)
        light1dot.move_to(light1position.get_value() * RIGHT)
        light2dot.move_to(light2position.get_value() * RIGHT)
        light1.move_to(light1dot.get_center())
        light2.move_to(light2dot.get_center())
        verticaldot_position = ValueTracker(5)
        verticaldot = Dot().move_to(verticaldot_position.get_value() * UP)
        self.play(
            FadeIn(light1),
            FadeIn(light2)
        )
        self.wait()
        light1dot.add_updater(lambda obj: obj.move_to(light1position.get_value() * RIGHT))
        light2dot.add_updater(lambda obj: obj.move_to(light2position.get_value() * RIGHT))
        light1.add_updater(lambda obj: obj.move_to(light1dot.get_center()))
        light2.add_updater(lambda obj: obj.move_to(light2dot.get_center()))
        self.add(light1, light2, light1dot, light2dot)
        self.play(Write(verticaldot))
        verticaldot.add_updater(lambda obj: obj.move_to(verticaldot_position.get_value() * UP))
        self.add(verticaldot)
        self.wait()
        line1 = Line(light1dot.get_center(), verticaldot.get_center(), color = BLUE)
        line2 = Line(light2dot.get_center(), verticaldot.get_center(), color = BLUE)
        self.play(Write(line1), Write(line2))
        self.wait()
        line1.add_updater(lambda obj: obj.become(Line(light1dot.get_center(), verticaldot.get_center(), color = BLUE)))
        line2.add_updater(lambda obj: obj.become(Line(light2dot.get_center(), verticaldot.get_center(), color = BLUE)))
        self.add(line1, line2)
        def get_perp_pos(a0, b0):
            a, b = a0, b0
            c2 = a * a + b * b
            px = b * a * b / c2
            py = b * a * a / c2
            return px * RIGHT + py * UP
        perp1 = Dot().set_color(YELLOW).move_to(get_perp_pos(light1position.get_value(), verticaldot_position.get_value()))
        perp2 = Dot().set_color(YELLOW).move_to(get_perp_pos(light2position.get_value(), verticaldot_position.get_value()))
        perpline1 = Line(ORIGIN, perp1.get_center(), color = GREEN)
        perpline2 = Line(ORIGIN, perp2.get_center(), color = GREEN)
        org = Dot()

        def right_angle(obj):
            obj_center = normalize(obj.get_center())
            a, b = obj_center[0], obj_center[1]
            res = VMobject(color = WHITE, fill_opacity = 1 / 2)
            res.set_points_as_corners(
                [
                    ORIGIN,
                    a * RIGHT + b * UP,
                    (a - b) * RIGHT + (a + b) * UP,
                    -b * RIGHT + a * UP,
                    ORIGIN
                ]
            )
            res.shift(obj.get_center())
            res.rotate(PI / 2, about_point = obj.get_center())
            res.scale(3 / 8, about_point = obj.get_center())
            return res
        
        right1, right2 = right_angle(perp1), right_angle(perp2)

        self.play(
            Write(org),
            Write(perp1),
            Write(perp2),
            Write(perpline1),
            Write(perpline2),
            Write(right1),
            Write(right2)
        )
        self.wait()

        perp1.add_updater(lambda obj: obj.move_to(get_perp_pos(light1position.get_value(), verticaldot_position.get_value())))
        perp2.add_updater(lambda obj: obj.move_to(get_perp_pos(light2position.get_value(), verticaldot_position.get_value())))
        perpline1.add_updater(lambda obj: obj.become(Line(ORIGIN, perp1.get_center(), color = GREEN)))
        perpline2.add_updater(lambda obj: obj.become(Line(ORIGIN, perp2.get_center(), color = GREEN)))
        right1.add_updater(lambda obj: obj.become(right_angle(perp1)))
        right2.add_updater(lambda obj: obj.become(right_angle(perp2)))
        self.add(perp1, perp2, perpline1, perpline2, right1, right2)
        def harmonic_position(a, b):
            if a == 0 or b == 0:
                return ORIGIN
            c = 1 / a + 1 / b
            c = 1 / c
            return c * RIGHT
        harmonic_dot = Dot(fill_opacity = 0)
        harmonic_dot.move_to(harmonic_position(light1position.get_value(), light2position.get_value()))
        harmonic_dot.add_updater(lambda obj: obj.move_to(harmonic_position(light1position.get_value(), light2position.get_value())))
        self.add(harmonic_dot)
        conn_line = VGroup(Line(perp1.get_center(), harmonic_dot.get_center()), Line(perp2.get_center(), harmonic_dot.get_center()))
        self.play(Write(conn_line))
        self.wait()
        conn_line.add_updater(lambda obj: obj.become(VGroup(Line(perp1.get_center(), harmonic_dot.get_center()), Line(perp2.get_center(), harmonic_dot.get_center()))))
        harmonic_light = lighthouse(col = BLUE).move_to(harmonic_dot.get_center())
        self.play(FadeIn(harmonic_light))
        harmonic_light.add_updater(lambda obj: obj.move_to(harmonic_dot.get_center()))
        self.play(light2position.animate.increment_value(1.5), run_time = 2.5)
        self.play(light1position.animate.increment_value(-1), run_time = 2)
        self.play(light1position.animate.increment_value(1.5), run_time = 2)
        self.play(verticaldot_position.animate.increment_value(-1.5), rate_func = there_and_back, run_time = 4)
        #self.play(light1position.animate.increment_value(-4), run_time = 2)
        self.wait()
        little_a, little_b, little_c = MathTex("a").set_color(RED), MathTex("b").set_color(GREEN), MathTex("c").set_color(YELLOW)
        for obj in [little_a, little_b, little_c]:
            obj.scale(1 + 1 / 4)
        little_a.next_to(light1dot, DOWN, buff = 1 / 2 + 1 / 4 + 1 / 8)
        little_b.next_to(light2dot, DOWN, buff = 1 / 2 + 1 / 4 + 1 / 8)
        little_c.next_to(harmonic_dot, DOWN, buff = 1 / 2 + 1 / 4 + 1 / 8)
        self.play(Write(little_a), Write(little_b), Write(little_c))
        little_a.add_updater(lambda obj: obj.next_to(light1dot, DOWN, buff = 1 / 2 + 1 / 4 + 1 / 8))
        little_b.add_updater(lambda obj: obj.next_to(light2dot, DOWN, buff = 1 / 2 + 1 / 4 + 1 / 8))
        little_c.add_updater(lambda obj: obj.next_to(harmonic_dot, DOWN, buff = 1 / 2 + 1 / 4 + 1 / 8))
        inverse_sum = VGroup(
            MathTex("\\frac{1}{a}"),
            MathTex("+\\frac{1}{b}"),
            MathTex("=\\frac{1}{c}")
        )
        for obj in inverse_sum:
            obj.scale(1 + 1 / 4)
        inverse_sum[0][0][-1].set_color(RED)
        inverse_sum[1][0][-1].set_color(GREEN)
        inverse_sum[2][0][-1].set_color(YELLOW)
        inverse_sum.arrange(RIGHT)
        inverse_sum.shift(8 * RIGHT + 4 * UP)
        self.play(TransformFromCopy(little_a, inverse_sum[0]))
        self.wait()
        self.play(TransformFromCopy(little_b, inverse_sum[1]))
        self.wait()
        self.play(TransformFromCopy(little_c, inverse_sum[2]))
        #self.wait()
        inverse_sum_surr = SurroundingRectangle(inverse_sum)
        self.play(Write(inverse_sum_surr))
        self.wait()
        self.play(
            light2position.animate.set_value(3),
            light1position.animate.set_value(-2),
            run_time = 4
        )
        self.play(
            self.camera.frame.animate.shift(6 * LEFT),
            VGroup(inverse_sum, inverse_sum_surr).animate.shift(13 * LEFT + DOWN),
        )
        self.wait()
        inv_sum_title = MathTex("\\text{Inverse Sum Theorem}").set_color(YELLOW).scale(3 / 2).shift(5.75 * UP)
        inv_sum_title.set_background_stroke(width = 5)
        self.play(Write(inv_sum_title))
        self.wait()
        three_proofs = VGroup(
            MathTex("\\text{Coordinate Geometry}"),
            MathTex("\\text{Trigonometry}"),
            MathTex("\\text{Circle Inversion}"),
        )
        blackscreen = Rectangle(width = 14, height = 10, fill_opacity = 0.75, fill_color = BLACK, stroke_width = 1 / 128)
        blackscreen.shift(LEFT + 2.5 * UP)
        three_proofs.scale(1.25)
        three_proofs[0].shift(4 * LEFT + 2.5 * UP)
        three_proofs[1].shift(-4 * LEFT - 2.5 * UP)
        three_proofs.shift(LEFT + 2.5 * UP)
        self.play(FadeIn(blackscreen))
        self.wait()
        self.play(Write(three_proofs[0], run_time = 2))
        self.wait()
        self.play(Write(three_proofs[1], run_time = 2))
        self.wait()
        self.play(Write(three_proofs[2], run_time = 2))
        self.wait()
        diff_levels = VGroup(
            MathTex("\\text{(straightforward)}"),
            MathTex("\\text{(easier to solve)}"),
            MathTex("\\text{(very elegant)}")
        ).set_color(YELLOW)
        diff_levels[0].next_to(three_proofs[0], RIGHT)
        diff_levels[1].next_to(three_proofs[1], LEFT)
        diff_levels[2].next_to(three_proofs[2], RIGHT)
        self.play(Write(diff_levels), run_time = 4.5)
        self.wait()
        self.play(
            FadeOut(diff_levels),
            FadeOut(three_proofs),
            FadeOut(blackscreen)
        )
        self.wait()
        twotoone = VGroup(
            MathTex("\\text{Two }", "\\text{lighthouses}"),
            MathTex("\\text{Single }", "\\text{lighthouse}"),
        )
        twotoone[0].shift(LEFT + 2.5 * UP + 4.5 * LEFT + UP)
        twotoone[0][1].next_to(twotoone[0][0], DOWN)
        twotoone[1][1].next_to(twotoone[1][0], DOWN)
        twotoone[1].next_to(twotoone[0], RIGHT)
        twotoone[1].shift(2 * RIGHT)
        for obj in twotoone:
            obj.set_background_stroke(width = 5)
        self.play(VGroup(inverse_sum, inverse_sum_surr).animate.shift(9 * RIGHT + UP))
        self.play(
            Write(twotoone),
            run_time = 3
        )
        twotoone_arrow = Arrow(twotoone[0].get_right(), twotoone[1].get_left())
        self.play(Write(twotoone_arrow))
        self.wait()
        self.play(twotoone_arrow.animate.rotate(540 * DEGREES), run_time = 2)
        self.wait()
        self.play(
            FadeOut(VGroup(inverse_sum, inverse_sum_surr), shift = RIGHT),
            FadeOut(VGroup(twotoone, twotoone_arrow), shift = LEFT),
            FadeOut(inv_sum_title, shift = UP),
            FadeOut(little_a), FadeOut(little_b), FadeOut(little_c)
        )
        self.wait()
        hidden_circle = Circle(radius = verticaldot_position.get_value() / 2).shift(verticaldot_position.get_value() / 2 * UP)
        self.play(Write(hidden_circle))
        self.wait()
        for obj in [line1, line2, perpline1, perpline2, conn_line, right1, right2, perp1, perp2, light1, light2]:
            obj.clear_updaters()
        self.play(
            FadeOut(VGroup(line1, line2, perpline1, perpline2, conn_line, right1, right2, perp1, perp2, light1, light2)),
        )
        self.wait()
        self.play(Write(conn_line))
        self.play(Write(VGroup(perp1, perp2)))
        self.wait()
        self.play(Write(VGroup(line1, line2)))
        self.play(FadeIn(VGroup(light1, light2)))
        self.wait()
        light1.add_updater(lambda obj: obj.move_to(light1dot.get_center()))
        light2.add_updater(lambda obj: obj.move_to(light2dot.get_center()))
        line1.add_updater(lambda obj: obj.become(Line(light1dot.get_center(), verticaldot.get_center(), color = BLUE)))
        line2.add_updater(lambda obj: obj.become(Line(light2dot.get_center(), verticaldot.get_center(), color = BLUE)))
        perp1.add_updater(lambda obj: obj.move_to(get_perp_pos(light1position.get_value(), verticaldot_position.get_value())))
        perp2.add_updater(lambda obj: obj.move_to(get_perp_pos(light2position.get_value(), verticaldot_position.get_value())))
        conn_line.add_updater(lambda obj: obj.become(VGroup(Line(perp1.get_center(), harmonic_dot.get_center()), Line(perp2.get_center(), harmonic_dot.get_center()))))
        self.add(light1, light2, line1, line2, perp1, perp2, conn_line)
        self.play(
            light1position.animate.increment_value(1),
            light2position.animate.increment_value(1),
            rate_func = there_and_back,
            run_time = 6
        )
        self.play(light2position.animate.increment_value(1.5), run_time = 3)
        self.wait()

class LakeAgain(MovingCameraScene):
    def construct(self):
        bub = SVGMobject("/Users/muthuveerappanr/Downloads/Manim/BaselProblem/lighthouse-svgrepo-com").set_color(WHITE)
        def lighthouse(rad = 10.0, col = YELLOW, fac = 1.0, lfac = 10):
            levels = int(lfac * rad)
            lights = VMobject()
            dr = rad / levels
            for r in np.arange(0, rad, dr):
                alpha = fac * 1.0 / (r + 1) ** 2
                ann = Annulus(
                    inner_radius = r,
                    outer_radius = r + dr,
                    color = col,
                    fill_opacity = alpha
                )
                lights.add(ann)
            lights.add(bub.copy().fade(3 / 8).scale(3 / 8).shift(DOWN / 2 + UP / 8))
            return lights
        lake = Circle(radius = 3.25, color = BLUE, fill_color = BLUE, fill_opacity = 1 / 8)
        observer = Dot().move_to(lake.get_bottom())
        self.play(Write(lake))
        self.play(Write(observer))
        self.wait()
        first_point = Dot(color = YELLOW)
        first_point.move_to(lake.get_bottom())
        first_point.rotate(2 * PI / 5, about_point = lake.get_center())
        diameter_line = Line(lake.get_top(), lake.get_bottom())
        top_connecting_line = Line(lake.get_top(), first_point.get_center())
        fivebypi = MathTex("\\frac{5}{\\pi}")
        dia_brace = Brace(diameter_line, LEFT)
        fivebypi.next_to(dia_brace, LEFT)
        pibyfive = MathTex("\\frac{\\pi}{5}")
        top_angle = Angle(diameter_line, top_connecting_line)
        pibyfive.next_to(top_angle, DOWN)
        pibyfive.shift(RIGHT / 4)
        pibyfive.set_background_stroke(width = 5)
        self.play(Write(first_point))
        self.wait()
        self.play(
            Write(diameter_line),
            Write(top_connecting_line)
        )
        self.play(
            Write(top_angle),
            Write(pibyfive)
        )
        self.wait()
        self.play(
            Write(dia_brace),
            Write(fivebypi)
        )
        self.wait()
        hor_line = Line(100 * LEFT, 10 * RIGHT, color = GREY).shift(3.25 * DOWN)
        top_connecting_line_extended = Line(lake.get_top(), 3.25 * DOWN + 2 * 3.25 * np.tan(PI / 5) * RIGHT)
        self.play(Write(hor_line), Write(top_connecting_line_extended))
        light1 = lighthouse(col = BLUE)
        light1.move_to(top_connecting_line_extended.get_end())
        self.play(FadeIn(light1))
        self.wait()
        self.play(
            FadeOut(pibyfive),
            FadeOut(fivebypi),
            FadeOut(dia_brace),
            FadeOut(diameter_line),
            FadeOut(top_angle),
            FadeOut(top_connecting_line),
            FadeOut(top_connecting_line_extended),
            FadeOut(first_point)
        )
        lake2 = Circle(radius = 3.25, color = BLUE, fill_color = BLUE, fill_opacity = 1 / 8)
        self.play(
            lake.animate.scale(1 / 2, about_point = 3.25 * DOWN),
            light1.animate.scale(1 / 2, about_point = 3.25 * DOWN),
            FadeIn(lake2)
        )
        self.wait()
        lake2_center = Dot(lake2.get_center())
        line2 = Line(3.25 * UP, 3.25 * DOWN / np.cos(PI / 5)).rotate_about_origin(PI / 5)
        ind21, ind22 = Dot(color = YELLOW), Dot(color = YELLOW)
        ind21.move_to(lake2.get_top())
        ind22.move_to(lake2.get_bottom())
        ind21.rotate(PI / 5, about_point = lake2.get_center())
        ind22.rotate(PI / 5, about_point = lake2.get_center())
        self.play(Write(lake2_center))
        self.wait()
        self.play(Write(line2))
        self.play(Write(VGroup(ind21, ind22)))
        self.wait()
        line21 = Line(lake2.get_top(), 3.25 * DOWN + 2 * 3.25 * np.tan(PI / 10) * RIGHT, color = GREEN)
        line22 = Line(lake2.get_top(), 3.25 * DOWN + 2 * 3.25 * np.tan(PI / 10 + PI / 2) * RIGHT, color = GREEN)
        self.play(
            Write(line21),
            Write(line22),
        )
        light21 = light1.copy().move_to(line21.get_end())
        light22 = light1.copy().move_to(line22.get_end())
        self.play(TransformFromCopy(light1, light21))
        self.play(
            TransformFromCopy(light1, light22),
            FadeOut(light1),
            self.camera.frame.animate.shift(15 * LEFT),
            run_time = 4
        )
        self.play(self.camera.frame.animate.shift(15 * RIGHT), run_time = 2)
        self.wait()
        self.play(
            FadeOut(line21),
            FadeOut(line22),
            FadeOut(line2),
            FadeOut(lake),
            FadeOut(ind21),
            FadeOut(ind22),
            FadeOut(lake2_center)
        )
        self.wait()
        light_group2 = VGroup()
        green_line_group2 = VGroup()
        ind_group2 = VGroup()
        line_group2 = VGroup()
        for k in range(4):
            light_group2.add(lighthouse(col = BLUE).scale(1 / 2).shift(3.25 * DOWN + 2 * 3.25 * np.tan((PI / 5 + k * PI) / 4) * RIGHT))
            green_line_group2.add(Line(3.25 * UP, light_group2[-1].get_center(), color = GREEN))
            ind_group2.add(Dot(color = YELLOW).move_to(3.25 * DOWN).rotate_about_origin(PI / 10 + k * PI / 2))
        
        self.play(
            lake2.animate.scale(1 / 2, about_point = 3.25 * DOWN),
            light21.animate.move_to(3.25 * DOWN + 2 * 3.25 * np.tan(PI / 10) * RIGHT / 2),
            light22.animate.move_to(3.25 * DOWN + 2 * 3.25 * np.tan(PI / 10 + PI / 2) * RIGHT / 2)
        )
        self.wait()
        line_group2.add(Line(light21.get_center(), ind_group2[2].get_center()))
        line_group2.add(Line(light22.get_center(), ind_group2[1].get_center()))
        lake3 = Circle(radius = 3.25, color = BLUE, fill_color = BLUE, fill_opacity = 1 / 8)
        self.play(Write(lake3))
        self.wait()
        #self.play(FadeIn(light_group2), Write(green_line_group2), Write(ind_group2), Write(line_group2))
        self.play(Write(line_group2), run_time = 2)
        self.wait()
        self.play(Write(ind_group2), run_time = 2)
        self.wait()
        self.play(Write(green_line_group2), run_time = 2)
        self.wait()
        self.play(FadeIn(light_group2), FadeOut(light21), FadeOut(light22))
        self.wait()
        self.play(
            FadeOut(line_group2),
            FadeOut(green_line_group2),
            FadeOut(lake2)
        )
        self.wait()
        def get_inds(gen):
            rad = 3.25 * 2 ** (gen - 2) / 2
            tcircle = Circle(radius = rad, color = BLUE, fill_color = BLUE, fill_opacity = 1 / 8)
            tcircle.move_to(3.25 * DOWN + tcircle.get_center() - tcircle.get_bottom())
            res = VGroup(tcircle)
            lights = VGroup()
            if gen > 4:
                mmin, mmax = -5, 5
            else:
                mmin, mmax = 0, 2 ** gen
            for k in range(mmin, mmax):
                res.add(Dot(color = YELLOW).move_to(tcircle.get_bottom()).rotate(2 * PI / 5 / 2 ** gen + k * 2 * PI / 2 ** gen, about_point = tcircle.get_center()))
                lights.add(lighthouse(col = BLUE).scale(1 / 2).move_to(3.25 * DOWN + 2 * rad * np.tan(PI / 5 / 2 ** gen + k * PI / 2 ** gen) * RIGHT))
            return VGroup(res, lights)
        ind3 = get_inds(3)
        self.play(
            FadeOut(VGroup(lake3, ind_group2)),
            FadeOut(light_group2),
            FadeIn(ind3),
        )
        self.wait()
        ind4 = get_inds(4)
        self.play(
            FadeOut(ind3),
            FadeIn(ind4),
        )
        self.wait()
        ind5 = get_inds(5)
        self.play(
            FadeOut(ind4),
            FadeIn(ind5),
        )
        self.wait()
        ind6 = get_inds(6)
        self.play(
            FadeOut(ind5),
            FadeIn(ind6),
        )
        self.wait()
        ind7 = get_inds(7)
        self.play(
            FadeOut(ind6),
            FadeIn(ind7),
        )
        self.wait()
        ind8 = get_inds(8)
        self.play(
            FadeOut(ind7),
            FadeIn(ind8),
        )
        self.wait()
        ind9 = get_inds(9)
        self.play(
            FadeOut(ind8),
            FadeIn(ind9),
        )
        self.wait()
        line_lighthouses = VGroup(hor_line.copy().set_color(BLUE))
        for k in range(-14, 16, 5):
            temp_lighthouse = lighthouse(col = BLUE).scale(1 / 2)
            #temp_lighthouse.move_to(observer.get_center())
            temp_lighthouse.shift(3.25 * PI / 20 * k * RIGHT)
            line_lighthouses.add(temp_lighthouse)
        line_lighthouses.shift(3.25 * DOWN)
        ticks = VGroup()
        tick = Line(UP / 8, DOWN / 8, color = GREEN)
        for k in range(-15, 23):
            j = k % 5
            if j == 1:
                continue
            ticks.add(tick.copy().shift(3.25 * DOWN + 3.25 * PI / 20 * k * RIGHT))
        self.play(
            FadeOut(ind9),
            FadeIn(line_lighthouses)
        )
        self.wait()
        self.play(ShowIncreasingSubsets(ticks), run_time = 3)
        self.wait()
        onecommazero = MathTex("(1,0)").set_background_stroke(width = 5).move_to(3.25 * DOWN + 3.25 * PI / 20 * RIGHT)
        onecommazero.shift(UR / 4)
        self.play(Write(onecommazero))
        self.wait()
        total_brightness = MathTex("\\cdots+\\frac{1}{-9}+\\frac{1}{-4}+\\frac{1}{1}+\\frac{1}{6}+\\frac{1}{11}+\\cdots")
        total_brightness.shift(2 * UP)
        total_brightness.set_background_stroke(width = 5)
        self.play(Write(total_brightness), run_time = 4)
        self.wait()
        original_lake = Circle(radius = 3.25, color = BLUE, fill_color = BLUE, fill_opacity = 1 / 8)
        self.play(
            total_brightness.animate.shift(4 * UP),
            FadeIn(pibyfive),
            FadeIn(fivebypi),
            FadeIn(dia_brace),
            FadeIn(diameter_line),
            FadeIn(top_angle),
            FadeIn(top_connecting_line),
            FadeIn(top_connecting_line_extended),
            FadeIn(first_point),
            FadeOut(onecommazero, shift = DOWN),
            FadeIn(original_lake),
            FadeOut(ticks),
            FadeOut(line_lighthouses)
        )
        total_line = Line(3.25 * DOWN, top_connecting_line_extended.get_end(), color = YELLOW)
        total_length = MathTex("\\frac{5}{\\pi}\\tan\\left(\\frac{\\pi}{5}\\right)")
        total_length.set_color(YELLOW)
        total_length.set_background_stroke(width = 5)
        total_length.shift(3.5 * DOWN + 2.125 * RIGHT + UP)
        self.play(Write(total_length), Write(total_line))
        self.wait()
        resulteqn = VGroup(
            MathTex("="),
            MathTex("1", "\\over", "\\frac{5}{\\pi}\\tan\\left(\\frac{\\pi}{5}\\right)")
        )
        resulteqn[1].next_to(resulteqn[0], LEFT)
        resulteqn.shift(2 * UP + 3 * LEFT)
        self.play(
            FadeOut(VGroup(total_line, pibyfive, fivebypi, dia_brace, diameter_line, top_angle, top_connecting_line, top_connecting_line_extended, first_point)),
            total_brightness.animate.next_to(resulteqn),
            Write(resulteqn[0]),
            Write(resulteqn[1][:-1]),
            total_length.animate.move_to(resulteqn[1][-1]).shift(DOWN / 4),
            FadeIn(line_lighthouses)
        )
        self.wait()
        resulteqn1 = MathTex("\\frac{\\pi}{\\tan\\left(\\frac{\\pi}{5}\\right)}", "=", "\\cdots+\\frac{1}{1/5-2}+\\frac{1}{1/5-1}+\\frac{1}{1/5}+\\frac{1}{1/5+1}+\\frac{1}{1/5+2}+\\cdots")
        resulteqn1.scale(1 - 1 / 8)
        resulteqn1.set_background_stroke(width = 5)
        resulteqn[1][-1].fade(1)
        #resulteqn1.shift(UP)
        self.play(Write(resulteqn1), run_time = 2)
        self.wait()
        self.play(
            VGroup(resulteqn, total_brightness, total_length).animate.shift(4 * UP),
            resulteqn1.animate.shift(2 * UP)
        )
        self.wait()
        resulteqn2 = MathTex("\\frac{\\pi}{\\tan (\\pi x)}=\\sum_{k \\in \\mathbb{Z}}\\frac{1}{x+k}")
        resulteqn2.set_background_stroke(width = 5)
        resulteqn2[0][7].set_color(YELLOW)
        resulteqn2[0][-3].set_color(YELLOW)
        resulteqn2[0][11].set_color(BLUE)
        resulteqn2[0][-1].set_color(BLUE)
        self.play(
            resulteqn1.animate.shift(4 * UP),
            Write(resulteqn2, run_time = 2),
            FadeOut(original_lake),
        )
        self.wait()

class Conclusion(MovingCameraScene):
    def construct(self):
        eqns = VGroup(
            MathTex("\\sin (\\pi x) = \\prod_{k \\in \\mathbb{Z}}(x+k)"),
            MathTex("\\frac{\\pi}{\\tan (\\pi x)}=\\sum_{k \\in \\mathbb{Z}}\\frac{1}{x+k}"),
            MathTex("\\frac{\\pi^2}{\\sin^2 (\\pi x)}=\\sum_{k \\in \\mathbb{Z}}\\frac{1}{(x+k)^2}"),
            MathTex("\\cdots"),
            MathTex("\\frac{(-1)^m}{m!}\\frac{d^m}{dx^m}\\left(\\frac{\\pi}{\\tan (\\pi x)}\\right)=\\sum_{k \\in \\mathbb{Z}}\\frac{1}{(x+k)^{m+1}}")
        )
        eqns.arrange(DOWN, buff = 1 - 1 / 4 - 1 / 8)
        eqns[-2].rotate(PI / 2)
        #eqns[0].set_color(BLUE)
        #eqns[1].set_color(YELLOW)
        #eqns[2].set_color(GREEN)
        #eqns[-1].set_color(TEAL)
        arrow0 = CurvedArrow(eqns[1].get_left() + LEFT / 2, eqns[0].get_left() + LEFT / 2, angle = -PI / 2)
        arrow1 = CurvedArrow(eqns[1].get_right() + RIGHT / 4, eqns[2].get_right() + RIGHT / 4, angle = -PI / 2)
        arrow2 = CurvedArrow(eqns[2].get_right() + RIGHT / 4, eqns[3].get_right() + 3 * RIGHT, angle = -PI / 2)
        arrow3 = CurvedArrow(eqns[3].get_right() + 3 * RIGHT, eqns[4].get_right() + RIGHT / 4, angle = -PI / 2)
        self.play(Write(eqns[1]), run_time = 2.5)
        self.wait()
        self.play(
            Write(eqns[2], run_time = 3),
            Write(arrow1),
        )
        self.wait()
        self.play(
            Write(eqns[3]),
            Write(eqns[-1], run_time = 3),
            #Write(arrow2),
            #Write(arrow3)
        )
        self.wait()
        self.play(
            Write(eqns[0], run_time = 2),
            Write(arrow0)
        )
        self.wait()
        bub = SVGMobject("/Users/muthuveerappanr/Downloads/Manim/BaselProblem/lighthouse-svgrepo-com").set_color(WHITE)
        def lighthouse(rad = 10.0, col = YELLOW, fac = 1.0, lfac = 10):
            levels = int(lfac * rad)
            lights = VMobject()
            dr = rad / levels
            for r in np.arange(0, rad, dr):
                alpha = fac * 1.0 / (r + 1) ** 2
                ann = Annulus(
                    inner_radius = r,
                    outer_radius = r + dr,
                    color = col,
                    fill_opacity = alpha
                )
                lights.add(ann)
            lights.add(bub.copy().fade(3 / 8).scale(3 / 8).shift(DOWN / 2 + UP / 8))
            return lights
        blue_val, yellow_val = ValueTracker(0), ValueTracker(PI / 8)
        def get_lights(start_val, col):
            res = VGroup()
            temp = Dot(fill_opacity = 0)
            temp.shift(5 * RIGHT)
            temp.rotate_about_origin(start_val)
            for k in range(8):
                res.add(lighthouse(col = col).scale(1 / 2).move_to(temp.get_center()))
                temp.rotate_about_origin(PI / 4)
            return res
        blue_lights = get_lights(blue_val.get_value(), col = BLUE)
        yellow_lights = get_lights(yellow_val.get_value(), col = YELLOW)
        self.play(
            AnimationGroup(FadeIn(blue_lights), FadeIn(yellow_lights)),
            run_time = 3
        )
        blue_lights.add_updater(lambda obj: obj.become(get_lights(blue_val.get_value(), col = BLUE)))
        yellow_lights.add_updater(lambda obj: obj.become(get_lights(yellow_val.get_value(), col = YELLOW)))
        self.add(blue_lights, yellow_lights)
        self.play(
            blue_val.animate.increment_value(3 * PI),
            yellow_val.animate.increment_value(-3 * PI),
            rate_func = linear,
            run_time = 10
        )
        self.wait()

if __name__ == "__main__":
    module_name = os.path.abspath(__file__)
    #output_location = "C:\ManimCE\media"
    #clear_cmd = "cls"
    #command_A = "manim " + module_name + " " + "RealCase" + " " + "-pql -n 42" + " --media_dir " + output_location
    #command_A = "manim " + module_name + " --media_dir " + output_location + " -pqh"
    #command_A = "manim "+ "-p" + " " + module_name + " " + "SphericalTautochrone" + " -n 0"+ " --renderer=opengl"
    command_A = "manim "+ "-pqh" + " " + module_name + " " + "OtherPointOnCircle" + " -n 0"
    #command_A = "manim "+ "-pqh" + " " + module_name + " " + "CircumferenceMethodProbability1" + " -n 0,5" + "  --disable_caching"
    #command_A = "manim "+ "-pqh" + " " + module_name + " --disable_caching"
    #command_A = "manim "+ "-pql" + " " + module_name
    #command_A = "manim "+ "-pqh" + " " + module_name
    #command_A = "manim "+ "--write_to_movie" + " " + module_name + " " + "SphericalTautochrone" + " -n 0"+ " --renderer=opengl"
    #os.system(clear_cmd)
    os.system(command_A)