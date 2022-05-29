from cmath import sqrt
from manim import *
from scipy import rand


class SquareToCircle(Scene):
    def construct(self):
        circle = Circle()  # create a circle
        circle.set_fill(PINK, opacity=0.5)  # set the color and transparency
        self.play(Create(circle))  # show the circle on screen

class HelloWorld(Scene):
    def construct(self):
        text = Text("Hello world", font_size=144)
        self.add(text)

class GlowingDot(Scene):
    def construct(self):
        dot = Dot(fill_opacity = 0.01)
        def get_glow(n, fac = 7):
            glow = VGroup()
            for k in range(n):
                temp = Circle(radius = (k + 1) / n, color = YELLOW, fill_opacity = (n - k) / n / fac, stroke_width = 0.01)
                glow.add(temp)
            return glow
        gl = get_glow(15)
        gl.add_updater(lambda z: z.move_to(dot))
        self.play(FadeIn(gl))
        self.add(gl)
        self.play(dot.animate.shift(UP))
        self.wait(0)

class MovingAround(Scene):
    def construct(self):
        n_lines = 100
        n = 20
        lines = VGroup()
        for k in range(n_lines):
            x = n * np.random.random()
            temp_line = Line(x * RIGHT + n * UP, x * RIGHT + n * DOWN)
            temp_line.rotate(2.0 * np.pi * np.random.random())
            temp_line.stroke_width = 0.75
            temp_line.set_opacity(0.75)
            lines.add(temp_line)
        self.play(ShowIncreasingSubsets(lines))
        self.wait()
        unit_circle = Circle()
        self.play(Write(unit_circle))

        def get_chord(line, circle):
            r = circle.get_radius()
            tangent = normalize(line.get_vector())
            normal = rotate_vector(tangent, PI / 2)
            center = circle.get_center()
            d1 = np.dot(line.get_start() - center, normal)
            if np.abs(d1) > r:
                return VectorizedPoint(center + r * normal)
            d2 = np.sqrt(r**2 - d1**2)
            chord = Line(
                center + d1 * normal - d2 * tangent,
                center + d1 * normal + d2 * tangent,
            )
            chord.color = GREEN
            chord.stroke_width = 0.75
            chord.stroke_opacity = 0.75
            return chord
        
        def create_chords():
            return VGroup(*[get_chord(line, unit_circle) for line in lines])
        
        chords = create_chords()
        self.play(FadeIn(chords))
        self.wait()

        chords.add_updater(lambda z: z.become(create_chords()))
        self.add(chords)
        self.play(unit_circle.animate.shift(LEFT))
        chords.clear_updaters()
        self.wait()

class TestGlow(Scene):
    def construct(self):
        rad = 3.5
        circ = Circle(radius = rad, color = WHITE)
        self.play(Write(circ))
        self.wait()

        chords = VGroup()
        n = 100
        for k in range(n):
            temp = rad * np.sqrt(np.random.random())
            chord_tip = np.sqrt(rad ** 2 - temp ** 2)
            temp_chord = Line(temp * RIGHT + chord_tip * UP, temp * RIGHT + chord_tip * DOWN)
            temp_chord.stroke_width = 0.5
            temp_chord.stroke_opacity = 0.5
            temp_chord.rotate_about_origin(np.random.random() * 2 * PI)
            chords.add(temp_chord)
        
        def get_glow(n, fac = 9):
            glow = VGroup()
            for k in range(n):
                temp = Circle(radius = (k + 1) / n, color = YELLOW, fill_opacity = (k + 1) / n / fac, stroke_width = 0.01)
                glow.add(temp)
            return glow
        gl = get_glow(15).scale(0.125)
        glows = VGroup(*[gl.copy().move_to(line.get_center()) for line in chords])
        self.wait()

        animlist = []
        for glo, ch in zip(glows, chords):
            animlist += [GrowFromCenter(glo), Write(ch), FadeOut(glo)]
        self.play(
            ShowIncreasingSubsets(chords),
            ShowSubmobjectsOneByOne(glows)
        )
        self.wait()

class RadialPoint(MovingCameraScene):
    def construct(self):
        rad = 3.5
        circ = Circle(radius = rad, color = WHITE)
        self.play(Write(circ))

        chords = VGroup()
        radial_lines = VGroup()
        n = 100
        for k in range(n):
            temp = rad * np.random.random()
            chord_tip = np.sqrt(rad ** 2 - temp ** 2)
            temp_chord = Line(temp * RIGHT + chord_tip * UP, temp * RIGHT + chord_tip * DOWN, color = BLUE)
            angle = np.random.random() * 2 * PI
            temp_radial_line = Line(ORIGIN, rad * RIGHT)
            temp_radial_line.stroke_width = 0.5
            temp_radial_line.stroke_color = YELLOW
            temp_dot = Dot(temp * RIGHT, color = YELLOW, stroke_width = 3, stroke_color = GREEN)
            radial_lines.add(VGroup(temp_chord, temp_radial_line, temp_dot).rotate_about_origin(angle))
            temp_chord_copy = temp_chord.copy()
            temp_chord_copy.stroke_width = 0.5
            temp_chord_copy.stroke_opacity = 0.5
            chords.add(temp_chord_copy)

        self.play(
            ShowIncreasingSubsets(chords),
            ShowSubmobjectsOneByOne(radial_lines),
            rate_func = linear,
            run_time = 10
        )
        self.play(FadeOut(radial_lines), run_time = 0.1)
        self.play(self.camera.frame.animate.shift(LEFT))
        self.wait()

class Introduction(MovingCameraScene):
    def construct(self):
        title = MathTex("\\text{Bertrand's Paradox}")
        title_line = Underline(title)
        title_group = VGroup(title, title_line)
        self.play(Write(title))
        self.play(Write(title_line))
        self.wait()
        self.play(title.animate.set_color(ORANGE))
        self.play(title_group.animate.to_edge(UP))
        self.wait()

        rad = 3
        circle = Circle(radius = rad, color = WHITE)
        inscribed_triangle = RegularPolygram(3, radius = rad, color = YELLOW)
        circle_group = VGroup(circle, inscribed_triangle).shift(0.5 * DOWN)
        self.play(Write(circle_group))
        self.wait()

        def get_chord_by_circumference(circle):
            circle_center = circle.get_center()
            circle_radius = circle.radius
            first_point = circle_center + circle_radius * DOWN
            second_point = Dot(first_point)
            second_point.rotate(2 * PI * np.random.random(), about_point = circle_center)
            second_point = second_point.get_center()
            chord = Line(first_point, second_point)
            chord.rotate(2 * PI * np.random.random(), about_point = circle_center)
            if chord.get_length() >= rad * np.sqrt(3):
                chord.color = BLUE
            else:
                chord.color = WHITE
            return chord
        
        def get_chord_by_radial_line(circle):
            circle_center = circle.get_center()
            circle_radius = circle.radius
            x = circle_radius * np.random.random()
            semi_chord = np.sqrt(circle_radius ** 2 - x * x)
            first_point = circle_center + x * RIGHT + semi_chord * UP
            second_point = circle_center + x * RIGHT + semi_chord * DOWN
            chord = Line(first_point, second_point)
            chord.rotate(2 * PI * np.random.random(), about_point = circle_center)
            if chord.get_length() >= rad * np.sqrt(3):
                chord.color = BLUE
            else:
                chord.color = WHITE
            return chord
        
        def get_chord_by_midpoint(circle):
            circle_center = circle.get_center()
            circle_radius = circle.radius
            r = circle_radius * np.sqrt(np.random.random())
            semichord = np.sqrt(circle_radius ** 2 - r * r)
            first_point = circle_center + r * RIGHT + semichord * UP
            second_point = circle_center + r * RIGHT + semichord * DOWN
            chord = Line(first_point, second_point)
            chord.rotate(2 * PI * np.random.random(), about_point = circle_center)
            if chord.get_length() >= rad * np.sqrt(3):
                chord.color = BLUE
            else:
                chord.color = WHITE
            return chord
        
        n = 200
        point_of_focus = 5 * DOWN
        second_circle_group = circle_group.copy().rotate(120 * DEGREES, about_point = point_of_focus)
        third_circle_group = circle_group.copy().rotate(-120 * DEGREES, about_point = point_of_focus)
        chords_by_circumference = VGroup(*[get_chord_by_circumference(circle_group[0]) for _ in range(n)])
        chords_by_radialline = VGroup(*[get_chord_by_radial_line(second_circle_group[0]) for _ in range(n)])
        chords_by_midpoint = VGroup(*[get_chord_by_midpoint(third_circle_group[0]) for _ in range(n)])
        chords_by_circumference_thicker = VGroup()
        chords_by_radialline_thicker = VGroup()
        chords_by_midpoint_thicker = VGroup()

        alternate_glow = VGroup(*[Line(ORIGIN, RIGHT / 4, color = YELLOW, stroke_width = 1 / 4).rotate_about_origin(k * 360 * DEGREES / 100) for k in range(100)])

        def add_end_dots(line):
            obj = line.copy()
            startingdot = Dot(obj.get_start(), color = YELLOW, stroke_color = RED, stroke_width = 3)
            endingdot = Dot(obj.get_end(), color = YELLOW, stroke_color = RED, stroke_width = 3)
            return VGroup(obj, startingdot, endingdot)
        
        for objc, objr, objm in zip(chords_by_circumference, chords_by_radialline, chords_by_midpoint):
            tobjc, tobjr, tobjm = objc.copy(), objr.copy(), objm.copy()
            for obj in [objc, objr, objm]:
                obj.stroke_width = 0.5
                obj.stroke_opacity = 0.35
            
            tobjc = add_end_dots(tobjc)
            tobjr = add_end_dots(tobjr)
            tobjm = add_end_dots(tobjm)

            chords_by_circumference_thicker.add(tobjc)

            temp_dot = Dot(tobjm[0].get_center(), color = LIGHT_PINK, stroke_color = YELLOW, stroke_width = 3)
            tobjm.add(alternate_glow.copy().move_to(tobjm[0].get_center()), temp_dot)
            chords_by_midpoint_thicker.add(tobjm)

            temp_radialline = Line(second_circle_group[0].get_center(), tobjr[0].get_center())
            temp_radialline = normalize(temp_radialline.get_vector())
            temp_radialline = Line(second_circle_group[0].get_center(), second_circle_group[0].get_center() + rad * temp_radialline, color = GREEN)
            tobjr.add(alternate_glow.copy().move_to(tobjr[0].get_center()), temp_radialline)
            chords_by_radialline_thicker.add(tobjr)
        
        self.add(second_circle_group, third_circle_group)
        self.play(
            self.camera.frame.animate.scale(7 / 4).move_to(point_of_focus + 1.75 * UP),
            FadeOut(title_group)
        )
        self.play(
            ShowIncreasingSubsets(chords_by_circumference),
            ShowSubmobjectsOneByOne(chords_by_circumference_thicker),
            ShowIncreasingSubsets(chords_by_radialline),
            ShowSubmobjectsOneByOne(chords_by_radialline_thicker),
            ShowIncreasingSubsets(chords_by_midpoint),
            ShowSubmobjectsOneByOne(chords_by_midpoint_thicker),
            rate_func = linear,
            run_time = 10
        )
        self.play(
            FadeOut(chords_by_circumference_thicker),
            FadeOut(chords_by_radialline_thicker),
            FadeOut(chords_by_midpoint_thicker),
            run_time = 0.01
        )

        self.wait()

class SimpleProblem(MovingCameraScene):
    def construct(self):
        title = MathTex("\\text{A Simpler Problem}")
        title_line = Underline(title)
        title_group = VGroup(title, title_line)
        self.play(Write(title))
        self.play(Write(title_line))
        self.play(title.animate.set_color(ORANGE))
        self.play(FadeOut(title_group, shift = UP))
        self.wait()

        rad = 3
        big_circle = Circle(radius = rad, color = WHITE)
        small_circle = Circle(radius = rad / 3, color = BLUE)
        self.play(
            GrowFromCenter(big_circle),
            #FadeIn(small_circle)
        )
        self.wait()

        def get_glow(n, fac = 7):
            glow = VGroup()
            for k in range(n):
                temp = Circle(radius = (k + 1) / n, color = YELLOW, fill_opacity = (n - k) / n / fac, stroke_width = 0.01)
                glow.add(temp)
            return glow

        n = 200
        glow = get_glow(15, fac = 9)
        dot_group, radius_group, glow_grp, random_group = VGroup(), VGroup(), VGroup(), VGroup()
        for _ in range(n):
            temp_group = VGroup()
            random_variable = np.sqrt(np.random.random()) * rad
            random_angle = 2 * PI * np.random.random()
            #temp_dot = Dot(random_variable * RIGHT, color = YELLOW, stroke_color = GREEN, stroke_width = 5)
            temp_dot = Dot(random_variable * RIGHT, radius = 1 / 32, color = GREEN)
            temp_dot.rotate_about_origin(random_angle)
            dot_group.add(temp_dot)
            random_radius = Line(ORIGIN, rad * RIGHT, stroke_color = random_color(), stroke_width = DEFAULT_STROKE_WIDTH / 2, stroke_opacity = 0.25)
            random_radius.rotate_about_origin(random_angle)
            radius_group.add(random_radius)
            random_glow = glow.copy().scale(0.375).move_to(temp_dot.get_center())
            glow_grp.add(random_glow)
            temp_group.add(random_radius)
            #temp_group.add(*random_glow)
            temp_group.add(temp_dot)
            random_group.add(temp_group)
        
        self.play(
            ShowIncreasingSubsets(dot_group),
            rate_func = linear,
            run_time = 10
        )
        self.wait()
        sub_radius_group = radius_group.copy()
        self.play(
            ShowSubmobjectsOneByOne(sub_radius_group),
            run_time = 6
        )
        self.play(FadeOut(sub_radius_group), run_time = 0.01)
        self.wait()

        self.play(LaggedStartMap(FadeIn, radius_group), run_time = 3)
        self.play(LaggedStartMap(FadeOut, radius_group), run_time = 3)
        self.wait()

        def check_points_in_circle(dot, circle):
            circle_radius = circle.radius
            distance_between = Line(dot.get_center(), circle.get_center()).get_length()
            if distance_between <= circle_radius:
                new_dot = Dot(dot.get_center(), color = BLUE)
                return new_dot
            else:
                return VectorizedPoint(ORIGIN)
        
        def points_inside_circle(circle):
            return(VGroup(*[check_points_in_circle(obj, circle) for obj in dot_group]))
        
        dots_inside = points_inside_circle(small_circle)
        self.play(
            FadeIn(dots_inside),
            FadeIn(small_circle)
        )
        #self.wait()

        dots_inside.add_updater(lambda z: z.become(points_inside_circle(small_circle)))
        self.add(dots_inside)
        self.play(small_circle.animate.shift(LEFT))
        self.play(Rotating(small_circle, radians = PI, about_point = ORIGIN, run_time = 3))
        self.play(small_circle.animate.move_to(UP + RIGHT), run_time = 1.5)
        dots_inside.clear_updaters()
        self.wait()

        statements = VGroup(
            MathTex("\\mathbb{P}", "(", "\\text{random point lying}\\\\", "\\text{inside a", " region}", ")"),
            MathTex("="),
            MathTex("\\text{Area}", "(", "\\text{region}", ")", "\\over", "\\text{Area}(", "\\text{Circle}", ")")
        )
        statements.arrange(DOWN, buff = 1)
        statements[0][4].set_color(BLUE)
        statements[-1][2].set_color(BLUE)
        statements[-1][-2].set_color(GREEN)
        statements.move_to(6 * RIGHT)

        self.play(
            self.camera.frame.animate.shift(3 * RIGHT),
            Write(statements[0])
        )
        self.wait()
        self.play(Write(statements[1:]))
        self.wait()

        probability_statements = VGroup(
            MathTex("\\mathbb{P}", "(", "\\text{random point lying inside}\\\\", "\\text{a circle of radius }", " r", "<1", ")"),
            MathTex("="),
            MathTex("\\pi \\cdot ", "r", "^2", "\\over", "\\pi \\cdot ", "1", "^2"),
            MathTex("r", "^2")
        )
        probability_statements[0].move_to(statements[0].get_center() - probability_statements[0].get_center())
        probability_statements[2].move_to(statements[2].get_center() - probability_statements[2].get_center())
        probability_statements[0][-3].set_color(BLUE)
        probability_statements[2][1].set_color(BLUE)
        probability_statements[2][-2].set_color(GREEN)

        strikeline = Line(ORIGIN, 0.75 * RIGHT, color = RED).rotate_about_origin(45 * DEGREES)
        strikeline_group = VGroup(
            strikeline.copy().move_to(probability_statements[2][0]),
            strikeline.copy().move_to(probability_statements[2][-3]),
        )
        strikeline_group.shift(LEFT / 16)

        self.play(
            FadeOut(dots_inside),
            small_circle.animate.move_to(ORIGIN).set_fill(BLUE),
            FadeTransform(statements[0], probability_statements[0]),
            FadeTransform(statements[2], probability_statements[2])
        )
        self.play(Write(strikeline_group))
        self.wait()

class CircumferenceMethod(MovingCameraScene):
    def construct(self):
        rad = 3.75
        pentagon_rad = 7
        base_circle = Circle(radius = rad, color = WHITE).shift(pentagon_rad * RIGHT)
        circle_group = VGroup()
        for k in range(5):
            circle_group.add(base_circle.copy().rotate_about_origin(k * 72 * DEGREES))
        #self.camera.frame.move_to(circle_group[0])
        self.camera.frame.scale(2)
        self.play(Write(circle_group))
        self.wait()

        method_names = VGroup(
            MathTex("\\text{Circumference method}"),
            MathTex("\\text{Radial line method}"),
            MathTex("\\text{Midpoint method}"),
            MathTex("\\text{Random point on circle,}\\\\", "\\text{line at random angle}"),
            MathTex("\\text{Random point on random radii,}\\\\", "\\text{line at random angle}"),
        )

        method_names[0].next_to(circle_group[0]).shift(RIGHT)
        method_names[1].next_to(circle_group[1]).shift(RIGHT)
        method_names[2].next_to(circle_group[2], direction = LEFT).shift(LEFT)
        method_names[3].next_to(circle_group[3], direction = LEFT).shift(LEFT)
        method_names[4].next_to(circle_group[4]).shift(RIGHT)

        self.play(
            self.camera.frame.animate.move_to(circle_group[0].get_center() + 3 * RIGHT).scale(1 / 2),
            Write(method_names)
        )
        self.play(method_names[0].animate.to_edge(UP).shift(0.5 * LEFT))
        self.wait()

        glow = VGroup(*[Line(ORIGIN, 1 * RIGHT / 4, color = YELLOW, stroke_width = 1 / 4).rotate_about_origin(k * 360 * DEGREES / 100) for k in range(100)])
        def get_chord_by_circumference(circle):
            circle_center = circle.get_center()
            circle_radius = circle.radius
            first_point = circle_center + circle_radius * DOWN
            second_point = Dot(first_point)
            second_point.rotate(2 * PI * np.random.random(), about_point = circle_center)
            second_point = second_point.get_center()
            chord = Line(first_point, second_point)
            chord.rotate(2 * PI * np.random.random(), about_point = circle_center)
            if chord.get_length() >= rad * np.sqrt(3):
                chord.color = BLUE
            else:
                chord.color = WHITE
            return chord
        
        def add_end_dots(line):
            obj = line.copy()
            startingdot = Dot(obj.get_start(), color = YELLOW, stroke_color = RED, stroke_width = 3)
            endingdot = Dot(obj.get_end(), color = YELLOW, stroke_color = RED, stroke_width = 3)
            return VGroup(obj, startingdot, endingdot)
        
        n = 1000
        chords = VGroup(*[get_chord_by_circumference(circle_group[0]) for _ in range(n)])
        thicker_chords = VGroup()
        for obj in chords:
            temp_group = obj.copy()
            obj.stroke_width = 0.3125
            obj.stroke_opacity = 0.25
            first_glow = glow.copy().move_to(obj.get_start())
            second_glow = glow.copy().move_to(obj.get_end())
            temp_group = add_end_dots(temp_group)
            temp_group.add(first_glow, second_glow)
            thicker_chords.add(temp_group)
        
        self.play(
            ShowIncreasingSubsets(chords),
            ShowSubmobjectsOneByOne(thicker_chords),
            rate_func = linear,
            run_time = 10
        )
        self.play(FadeOut(thicker_chords), run_time = 0.1)
        self.wait()

        new_chords = chords.copy()
        animlist = []
        for obj in chords:
            angle = Line(circle_group[0].get_center(), obj.get_start()).get_angle()
            animlist += [Rotate(obj, -angle - PI / 2, about_point = circle_group[0].get_center())]
        self.play(LaggedStart(*animlist), rate_func = linear, run_time = 5)
        self.wait()

        center_line = Line(circle_group[0].get_bottom() + 0.25 * DOWN, circle_group[0].get_top() + 0.25 * UP, color = BLUE)
        fliplist = []
        for obj in chords:
            if obj.get_start()[0] > obj.get_end()[0]:
                fliplist += [Rotate(obj, PI, axis = UP, about_point = center_line.get_start())]
        self.play(Write(center_line))
        self.play(LaggedStart(*fliplist), run_time = 5)
        self.wait()

        random_chord = Dot(circle_group[0].get_bottom())
        random_chord.rotate(PI * 3 / 4, about_point = circle_group[0].get_center())
        random_chord = Line(circle_group[0].get_bottom(), random_chord.get_center(), color = ORANGE)
        random_chord = add_end_dots(random_chord)

        tangent_line = Line(circle_group[0].get_bottom() + LEFT, circle_group[0].get_bottom() + RIGHT)
        tangent_angle = Angle(tangent_line, random_chord[0], radius = 0.5)
        theta_text = MathTex("\\theta").set_background_stroke(width = 5)
        theta_text.next_to(tangent_angle)
        self.play(Write(tangent_line), Write(random_chord))
        self.play(Write(tangent_angle), Write(theta_text))
        self.wait()

        theta_distribution = MathTex("\\Theta \\sim U(0,\\pi/2)")
        theta_distribution.move_to(circle_group[0].get_center())
        theta_distribution.shift(2 * UP + 7 * RIGHT)
        self.play(Write(theta_distribution))
        self.wait()

        center_to_end_line = Line(circle_group[0].get_center(), random_chord[0].get_end())
        center_to_end_line.stroke_width = 0.5
        center_angle = Angle(Line(center_line.get_center(), center_line.get_start()), center_to_end_line, radius = 0.5, other_angle = False)
        twotheta_text = MathTex("2\\theta").set_background_stroke(width = 5)
        twotheta_text.next_to(center_angle)
        center_dot = Dot(circle_group[0].get_center(), color = GREEN, stroke_color = BLUE, stroke_width = 3)
        self.play(Write(center_dot), Write(center_to_end_line), Write(center_angle), Write(twotheta_text))
        self.wait()

        chord_length = MathTex("L=2\\sin\\theta").set_background_stroke(width = 0.5)
        chord_length.move_to(random_chord.get_center() + 0.25 * LEFT)
        chord_length.rotate(PI * 3 / 4)

        expected_length = VGroup(
            MathTex("\\mathbb{E}(L)", "=", "\\int_0^{\\pi/2}2\\sin\\theta\\,\\frac{d\\theta}{\\pi/2}"),
            MathTex("=", "\\frac{4}{\\pi}"),
            MathTex("\\approx", "1.2732")
        )
        expected_length[0].move_to(circle_group[0].get_center() + 1 * UP / 4 + 7 * RIGHT)
        expected_length[1].move_to(expected_length[0][1].get_center() + expected_length[1].get_center() - expected_length[1][0].get_center())
        expected_length[2].move_to(expected_length[0][1].get_center() + expected_length[2].get_center() - expected_length[2][0].get_center())
        #expected_length.arrange(DOWN, buff = MED_LARGE_BUFF)
        expected_length[1].shift(3 * DOWN / 2)
        expected_length[2].shift(6 * DOWN / 2)
        self.play(Write(expected_length[0]))
        self.wait()
        self.play(Write(expected_length[1:]))
        self.wait()

        self.play(
            self.camera.frame.animate.move_to(circle_group[1].get_center() + 3 * RIGHT),
            FadeOut(method_names[0]),
            FadeOut(center_line),
            FadeOut(chords),
            FadeOut(random_chord)
        )
        self.wait()

class RadiallineMethod(MovingCameraScene):
    def construct(self):
        rad = 3.75
        pentagon_rad = 7
        base_circle = Circle(radius = rad, color = WHITE).shift(pentagon_rad * RIGHT)
        circle_group = VGroup()
        for k in range(5):
            circle_group.add(base_circle.copy().rotate_about_origin(k * 72 * DEGREES))
        self.camera.frame.move_to(circle_group[1].get_center() + 3 * RIGHT)

        method_names = VGroup(
            MathTex("\\text{Circumference method}"),
            MathTex("\\text{Radial line method}"),
            MathTex("\\text{Midpoint method}"),
            MathTex("\\text{Random point on circle,}\\\\", "\\text{line at random angle}"),
            MathTex("\\text{Random point on random radii,}\\\\", "\\text{line at random angle}"),
        )

        method_names[0].next_to(circle_group[0]).shift(RIGHT)
        method_names[1].next_to(circle_group[1]).shift(RIGHT)
        method_names[2].next_to(circle_group[2], direction = LEFT).shift(LEFT)
        method_names[3].next_to(circle_group[3], direction = LEFT).shift(LEFT)
        method_names[4].next_to(circle_group[4]).shift(RIGHT)

        self.add(circle_group, method_names)
        self.play(method_names[1].animate.shift(3.5 * UP))
        self.wait()

        n = 1000
        def get_chord_by_radial_line(circle):
            circle_center = circle.get_center()
            circle_radius = circle.radius
            x = circle_radius * np.random.random()
            semi_chord = np.sqrt(circle_radius ** 2 - x * x)
            first_point = circle_center + x * RIGHT + semi_chord * UP
            second_point = circle_center + x * RIGHT + semi_chord * DOWN
            chord = Line(first_point, second_point)
            chord.rotate(2 * PI * np.random.random(), about_point = circle_center)
            if chord.get_length() >= rad * np.sqrt(3):
                chord.color = BLUE
            else:
                chord.color = WHITE
            return chord
        
        glow = VGroup(*[Line(ORIGIN, 1 * RIGHT / 4, color = YELLOW, stroke_width = 1 / 4).rotate_about_origin(k * 360 * DEGREES / 100) for k in range(100)])
        chords = VGroup(*[get_chord_by_radial_line(circle_group[1]) for _ in range(n)])
        
        thicker_chords = VGroup()
        for obj in chords:
            temp_group = VGroup()
            temp_group.add(obj.copy())
            obj.stroke_width = 0.3125
            obj.stroke_opacity = 0.25
            random_radius = rad * normalize(Line(circle_group[1].get_center(), obj.get_center()).get_vector())
            random_radius = Line(circle_group[1].get_center(), circle_group[1].get_center() + random_radius, color = GREEN)
            first_glow = glow.copy().move_to(obj.get_center())
            temp_group.add(random_radius, first_glow)
            thicker_chords.add(temp_group)
        
        self.play(
            ShowIncreasingSubsets(chords),
            ShowSubmobjectsOneByOne(thicker_chords),
            run_time = 10,
            rate_func = linear
        )
        self.play(FadeOut(thicker_chords), run_time = 0.1)
        self.wait()

        horizontal_radius = Line(circle_group[1].get_center(), circle_group[1].get_center() + rad * RIGHT, color = GREEN)
        x = rad / 2
        horizontal_glow = glow.copy().shift(circle_group[1].get_center() + x * RIGHT)
        y = np.sqrt(rad ** 2 - x ** 2)
        vertical_chord = Line(x * RIGHT + y * UP, x * RIGHT + y * DOWN, color = BLUE).shift(circle_group[1].get_center())
        horizontal_brace = Brace(Line(ORIGIN, x * RIGHT)).shift(circle_group[1].get_center())
        smallx = MathTex("x").set_background_stroke(width = 5)
        center_dot = Dot(circle_group[1].get_center(), color = GREEN, stroke_color = BLUE, stroke_width = 3)
        smallx.next_to(horizontal_brace, DOWN)
        chord_length = MathTex("2\\sqrt{1-x^2}")
        chord_length.move_to(circle_group[1].get_center() + x * RIGHT).rotate(- PI / 2)
        chord_length.shift(RIGHT / 2)
        chord_length.set_background_stroke(width = 5)

        self.play(
            Write(horizontal_radius),
            Write(horizontal_glow),
            FadeIn(vertical_chord),
            Write(center_dot)
        )
        self.wait()
        self.play(
            Write(horizontal_brace),
            Write(smallx)
        )
        self.play(Write(chord_length))
        self.wait()

        expected_length = VGroup(
            MathTex("\\mathbb{E}(L)", "=", "\\int_0^1 2\\sqrt{1-x^2}\\,dx"),
            MathTex("=", "\\pi/2 \\approx 1.5708")
        )
        expected_length.move_to(circle_group[1].get_center() + 8 * UP / 4 + 7 * RIGHT)
        for obj in expected_length:
            obj.scale(0.875)
        expected_length[1].move_to(expected_length[0][1].get_center() + expected_length[1].get_center() - expected_length[1][0].get_center())
        expected_length[1].shift(DOWN)

        self.play(Write(expected_length[0]))
        self.wait()
        self.play(Write(expected_length[1]))
        self.wait()

        self.play(
            FadeOut(VGroup(horizontal_brace, horizontal_glow, horizontal_radius, smallx, chord_length, center_dot)),
            FadeOut(method_names[1], shift = UP),
            expected_length.animate.shift(UP)
        )
        self.wait()

        substitution = MathTex("x=\\cos\\theta").scale(0.875).move_to(expected_length[0]).shift(3 * DOWN)
        new_integral = MathTex("\\mathbb{E}(L)", "=", "\\int_0^{\\pi/2}2\\sin\\theta\\cdot", "\\sin\\theta", "\\,d\\theta")
        new_integral.scale(0.875).move_to(expected_length[0]).shift(4 * DOWN)

        self.play(Write(substitution))
        self.play(Write(new_integral))
        self.wait()

        animlist = []
        for obj in chords:
            angle = Line(circle_group[1].get_center(), obj.get_end()).get_angle()
            animlist += [Rotate(obj, -angle - PI / 2, about_point = circle_group[1].get_center())]
        angle = Line(circle_group[1].get_center(), vertical_chord.get_end()).get_angle()
        self.play(Rotate(vertical_chord, -angle - PI / 2, about_point = circle_group[1].get_center()))
        self.play(LaggedStart(*animlist), rate_func = linear, run_time = 5)
        self.wait()

        tangent_line = Line(vertical_chord.get_end() + 2 * LEFT, vertical_chord.get_end() + 2 * RIGHT, color = GREEN)
        center_line = Line(circle_group[1].get_center() + rad * DOWN + DOWN / 4, circle_group[1].get_center() + rad * UP + UP / 4, color = RED)
        tangent_angle = Angle(tangent_line, Line(vertical_chord.get_end(), vertical_chord.get_start()), radius = 0.75)
        theta_text = MathTex("\\theta").set_background_stroke(width = 5)
        theta_text.next_to(tangent_angle)
        theta_text.shift(UP / 4)
        twotheta_line = Line(circle_group[1].get_center(), vertical_chord.get_start(), stroke_width = 0.5)
        twotheta_angle = Angle(Line(circle_group[1].get_center(), center_line.get_start()), Line(circle_group[1].get_center(), vertical_chord.get_start()), radius = 0.5)
        twotheta_text = MathTex("2\\theta").next_to(twotheta_angle)
        twotheta_text.set_background_stroke(width = 5)
        self.play(Write(VGroup(tangent_line, center_line, center_dot, twotheta_line, twotheta_angle, twotheta_text, tangent_angle, theta_text)), run_time = 2)
        self.wait()

        theta_distribution = MathTex("f_{\\Theta}(\\theta)", "=", "\\sin\\theta")
        theta_distribution.move_to(expected_length)
        theta_distribution.shift(2 * DOWN)
        self.play(
            FadeOut(expected_length),
            FadeOut(substitution),
            new_integral.animate.move_to(expected_length)
        )
        box = SurroundingRectangle(new_integral[-2])
        self.wait()
        self.play(Write(theta_distribution), Create(box))
        self.wait()
        self.play(
            FadeOut(VGroup(vertical_chord, chords, tangent_line, center_line, center_dot, twotheta_line, twotheta_angle, twotheta_text, tangent_angle, theta_text)),
            self.camera.frame.animate.move_to(circle_group[2].get_center() + 3 * LEFT)
        )
        self.wait()

class ExpectedValue(MovingCameraScene):
    def construct(self):
        expected_value = VGroup(
            MathTex("\\text{expected value of }\\\\", "\\text{a }", "\\text{random }", "\\text{quantity}"),
            Tex("="),
            MathTex("\\int"),
            Tex("that ", "quantity"),
            MathTex("\\times"),
            MathTex("\\text{a measure of the}\\\\", "\\text{randomness}")
        )
        expected_value.arrange()
        expected_value[-1][-1].shift(3 * LEFT / 4)
        #expected_value[0][-2].set_color(YELLOW)
        #expected_value[0][-1].set_color(RED)
        expected_value[3][-1].set_color(RED)
        expected_value[-1][-1].set_color(YELLOW)
        self.play(Write(expected_value[0]))
        self.play(
            Write(expected_value[1]),
            expected_value[0][-2].animate.set_color(YELLOW),
            expected_value[0][-1].animate.set_color(RED)
        )
        self.wait()
        self.play(Write(expected_value[3]))
        self.play(Write(expected_value[4]))
        self.play(Write(expected_value[5]))
        self.play(Write(expected_value[2]))
        self.wait()

        box = SurroundingRectangle(expected_value[-1])
        density_function = Tex("Probability Density Function").shift(2 * UP)
        connecting_line = Line(box.get_top(), density_function.get_bottom())

        self.play(Write(box))
        self.wait()
        self.play(
            Write(density_function),
            Write(connecting_line)
        )
        self.wait()
        self.play(FadeOut(expected_value), FadeOut(density_function), FadeOut(connecting_line), FadeOut(box))
        self.wait()

        table = VGroup(
            MathTex("\\text{Random method}"), MathTex("\\text{Density function}"), MathTex("\\text{Expected value}"),
            MathTex("\\text{Circumference}"), MathTex("\\frac{2}{\\pi}"), MathTex("\\frac{4}{\\pi}"),
            MathTex("\\text{Radial line}"), MathTex("\\sin\\theta"), MathTex("\\frac{\\pi}{2}"),
            MathTex("\\text{Midpoint}"), MathTex("\\sin2\\theta"), MathTex("\\frac{4}{3}"),
            MathTex("\\text{Point on circle,} \\\\ \\text{line at random angle}"), MathTex("\\frac{4}{\\pi}\\sin^2\\theta"), MathTex("\\frac{16}{3\\pi}"),
            MathTex("\\text{Point on radial line,} \\\\ \\text{line at random angle}"), MathTex("\\frac{2}{\\pi}\\sin\\theta\\sinh^{-1}(\\tan\\theta)"), MathTex("\\frac{2}{\\pi}(2G+1)")
        )
        table.arrange_in_grid(6, 3)
        for k in [0, 3, 6, 9, 12, 15]:
            table[k].shift(LEFT / 2)
        for k in [0, 3, 6, 9, 12, 15]:
            table[k + 2].shift(RIGHT / 2)
        row_separators = VGroup(*[Line(12 * LEFT, 1.5 * RIGHT, color = BLUE) for _ in range(5)])
        for k in range(4):
            row_separators[k + 1].shift((table[3 * k + 5].get_bottom() + table[3 * k + 8].get_top()) / 2)
        row_separators[0].set_color(YELLOW)
        row_separators[0].next_to(table[1], direction = DOWN, buff = 1 / 8)
        row_separators[0].shift(LEFT / 2)
        temp_separator = row_separators[0].copy()
        temp_separator.shift(1 / 16)

        self.play(
            Write(table[:3]),
            Write(temp_separator), Write(row_separators[0]),
            LaggedStartMap(Create, row_separators[1:], run_time = 2)
        )
        for k in range(1, 6):
            self.play(Write(table[3 * k: 3 * k + 3]))
        self.wait()
        self.play(
            FadeOut(VGroup(table))
        )
        self.wait()
        table = VGroup(
            MathTex("\\text{Random method}"), MathTex("\\mathbb{E}(\\text{length})"), MathTex("\\mathbb{P}(\\text{pts. on opp. sides})"),
            MathTex("\\text{Circumference}"), MathTex("\\frac{4}{\\pi}"), MathTex("\\frac{1}{3}-\\frac{5}{4\\pi^2}"),
            MathTex("\\text{Radial line}"), MathTex("\\frac{\\pi}{2}"), MathTex("\\frac{128}{45\\pi^2}"),
            MathTex("\\text{Midpoint}"), MathTex("\\frac{4}{3}"), MathTex("\\frac{1}{8}+\\frac{2}{3\\pi^2}"),
            MathTex("\\text{Point on circle,} \\\\ \\text{line at random angle}"), MathTex("\\frac{16}{3\\pi}"), MathTex("\\frac{1}{3}"),
            MathTex("\\text{Point on radial line,} \\\\ \\text{line at random angle}"), MathTex("\\frac{2}{\\pi}(2G+1)"), MathTex("\\frac{427+60\\pi^2-480\\ln 2}{180\pi^2}")
        )
        table.arrange_in_grid(6, 3)
        for k in [0, 3, 6, 9, 12, 15]:
            table[k].shift(LEFT / 2)
        for k in [0, 3, 6, 9, 12, 15]:
            table[k + 2].shift(RIGHT / 2)
        for k in range(6):
            self.play(Write(table[3 * k: 3 * k + 3]), run_time = 3)
        self.wait()

class MidpointMethod(MovingCameraScene):
    def construct(self):
        rad = 3.75
        pentagon_rad = 7
        base_circle = Circle(radius = rad, color = WHITE).shift(pentagon_rad * RIGHT)
        circle_group = VGroup()
        for k in range(5):
            circle_group.add(base_circle.copy().rotate_about_origin(k * 72 * DEGREES))
        self.camera.frame.move_to(circle_group[2].get_center() + 3 * LEFT)

        method_names = VGroup(
            MathTex("\\text{Circumference method}"),
            MathTex("\\text{Radial line method}"),
            MathTex("\\text{Midpoint method}"),
            MathTex("\\text{Random point on circle,}\\\\", "\\text{line at random angle}"),
            MathTex("\\text{Random point on random radii,}\\\\", "\\text{line at random angle}"),
        )

        method_names[0].next_to(circle_group[0]).shift(RIGHT)
        method_names[1].next_to(circle_group[1]).shift(RIGHT)
        method_names[2].next_to(circle_group[2], direction = LEFT).shift(LEFT)
        method_names[3].next_to(circle_group[3], direction = LEFT).shift(LEFT)
        method_names[4].next_to(circle_group[4]).shift(RIGHT)

        self.add(circle_group, method_names)
        self.play(method_names[2].animate.shift(3.5 * UP))
        self.wait()

        def get_chord_by_midpoint(circle):
            circle_center = circle.get_center()
            circle_radius = circle.radius
            r = circle_radius * np.sqrt(np.random.random())
            semichord = np.sqrt(circle_radius ** 2 - r * r)
            first_point = circle_center + r * RIGHT + semichord * UP
            second_point = circle_center + r * RIGHT + semichord * DOWN
            chord = Line(first_point, second_point).rotate(2 * PI * np.random.random(), about_point = circle_group[2].get_center())
            if chord.get_length() >= rad * np.sqrt(3):
                chord.color = BLUE
            else:
                chord.color = WHITE
            return chord
        
        n = 1000
        glow = VGroup(*[Line(ORIGIN, 1 * RIGHT / 4, color = YELLOW, stroke_width = 1 / 4).rotate_about_origin(k * 360 * DEGREES / 100) for k in range(100)])
        chords = VGroup(*[get_chord_by_midpoint(circle_group[2]) for _ in range(n)])
        thicker_chords = VGroup()
        for obj in chords:
            temp_group = VGroup()
            #obj.stroke_color = BLUE
            temp_group.add(obj.copy())
            obj.stroke_width = 0.3125
            obj.stroke_opacity = 0.25
            first_glow = glow.copy().move_to(obj.get_center())
            temp_group.add(first_glow)
            thicker_chords.add(temp_group)
        
        self.play(
            ShowIncreasingSubsets(chords),
            ShowSubmobjectsOneByOne(thicker_chords),
            run_time = 10,
            rate_func = linear
        )
        self.play(FadeOut(thicker_chords), run_time = 0.1)
        self.wait()

        horizontal_radius = Line(circle_group[2].get_center(), circle_group[2].get_center() + rad * RIGHT, color = GREEN, stroke_color = GREEN)
        x = rad / 2
        horizontal_glow = glow.copy().shift(circle_group[2].get_center() + x * RIGHT)
        y = np.sqrt(rad ** 2 - x ** 2)
        vertical_chord = Line(x * RIGHT + y * UP, x * RIGHT + y * DOWN, color = BLUE).shift(circle_group[2].get_center())
        horizontal_brace = Brace(Line(ORIGIN, x * RIGHT)).shift(circle_group[2].get_center())
        smallx = MathTex("r").set_background_stroke(width = 5)
        center_dot = Dot(circle_group[2].get_center(), color = GREEN, stroke_color = BLUE, stroke_width = 3)
        smallx.next_to(horizontal_brace, DOWN)
        chord_length = MathTex("2\\sqrt{1-r^2}")
        chord_length.move_to(circle_group[2].get_center() + x * RIGHT).rotate(- PI / 2)
        chord_length.shift(RIGHT / 2)
        chord_length.set_background_stroke(width = 5)

        self.play(
            Write(horizontal_radius),
            FadeIn(horizontal_glow),
            FadeIn(vertical_chord),
            Write(center_dot)
        )
        self.wait()
        self.play(
            Write(horizontal_brace),
            Write(smallx)
        )
        self.play(Write(chord_length))
        self.wait()

        density = VGroup(
            MathTex("r \\sim U(0,1)"),
            MathTex("F(r)", "=", "\\mathbb{P}(\\text{distance} \\leq r)=r^2"),
            MathTex("f(r)", "=", "F'(r)=2r")
        )
        for obj in density:
            obj.scale(0.875)
        density.arrange(direction = DOWN)
        density.move_to(circle_group[2].get_center() + 8 * UP / 4 + 7 * LEFT)

        strikeline = Cross(density[0])

        self.play(Write(density[0]))
        self.wait()
        self.play(Write(strikeline))
        self.wait()
        self.play(Write(density[1]))
        self.wait()
        self.play(Write(density[2]))
        self.wait()

        self.play(
            FadeOut(VGroup(density[:2], strikeline)),
            density[2].animate.move_to(density[0])
        )
        self.wait()

        expected_length = VGroup(
            MathTex("\\mathbb{E}(L)", "=", "\\int_0^1 2\\sqrt{1-r^2}\\cdot 2r\\,dr"),
            MathTex("=", "4/3 \\approx 1.3334")
        )
        expected_length.move_to(circle_group[2].get_center() + 6 * UP / 4 + 7 * LEFT)
        for obj in expected_length:
            obj.scale(0.875)
        expected_length[1].move_to(expected_length[0][1].get_center() + expected_length[1].get_center() - expected_length[1][0].get_center())
        expected_length[1].shift(DOWN)

        self.play(Write(expected_length[0]))
        self.wait()
        self.play(Write(expected_length[1]))
        self.wait()

        self.play(
            FadeOut(VGroup(horizontal_brace, horizontal_glow, horizontal_radius, smallx, chord_length, center_dot)),
            FadeOut(method_names[2], shift = UP),
            FadeOut(density[-1], shift = UP),
            expected_length.animate.shift(UP)
        )
        self.wait()

        substitution = MathTex("r=\\cos\\theta").scale(0.875).move_to(expected_length[0]).shift(3 * DOWN)
        new_integral = MathTex("\\mathbb{E}(L)", "=", "\\int_0^{\\pi/2}2\\sin\\theta\\cdot", "\\sin 2\\theta", "\\,d\\theta")
        new_integral.scale(0.875).move_to(expected_length[0]).shift(4 * DOWN)
        self.play(Write(substitution))
        self.play(Write(new_integral))
        self.wait()

        animlist = []
        for obj in chords:
            angle = Line(circle_group[2].get_center(), obj.get_end()).get_angle()
            animlist += [Rotate(obj, -angle - PI / 2, about_point = circle_group[2].get_center())]
        angle = Line(circle_group[2].get_center(), vertical_chord.get_end()).get_angle()

        self.play(Rotate(vertical_chord, -angle - PI / 2, about_point = circle_group[2].get_center()))
        self.play(LaggedStart(*animlist), rate_func = linear, run_time = 5)
        self.wait()

        tangent_line = Line(vertical_chord.get_end() + 2 * LEFT, vertical_chord.get_end() + 2 * RIGHT, color = GREEN)
        center_line = Line(circle_group[2].get_center() + rad * DOWN + DOWN / 4, circle_group[2].get_center() + rad * UP + UP / 4, color = RED)
        tangent_angle = Angle(tangent_line, Line(vertical_chord.get_end(), vertical_chord.get_start()), radius = 0.75)
        theta_text = MathTex("\\theta").set_background_stroke(width = 5)
        theta_text.next_to(tangent_angle)
        theta_text.shift(UP / 4)
        theta_distribution = MathTex("f_{\\Theta}(\\theta)", "=", "\\sin2\\theta")
        theta_distribution.move_to(new_integral)
        theta_distribution.shift(1.5 * DOWN)
        box = SurroundingRectangle(new_integral[-2])

        self.play(Write(VGroup(tangent_line, center_line, center_dot, tangent_angle, theta_text, theta_distribution, box)), run_time = 2)
        self.wait()
        self.play(
            FadeOut(VGroup(tangent_line, center_line, center_dot, tangent_angle, theta_text, chords, vertical_chord)),
            self.camera.frame.animate.move_to(circle_group[3].get_center() + 3 * LEFT)
        )
        self.wait()

class RandomPointRandomAngle3(MovingCameraScene):
    def construct(self):
        rad = 3.75
        pentagon_rad = 7
        base_circle = Circle(radius = rad, color = WHITE).shift(pentagon_rad * RIGHT)
        circle_group = VGroup()
        for k in range(5):
            circle_group.add(base_circle.copy().rotate_about_origin(k * 72 * DEGREES))
        self.camera.frame.move_to(circle_group[3].get_center() + 3 * LEFT)

        method_names = VGroup(
            MathTex("\\text{Circumference method}"),
            MathTex("\\text{Radial line method}"),
            MathTex("\\text{Midpoint method}"),
            MathTex("\\text{Random point on circle,}\\\\", "\\text{line at random angle}"),
            MathTex("\\text{Random point on random radii,}\\\\", "\\text{line at random angle}"),
        )

        method_names[0].next_to(circle_group[0]).shift(RIGHT)
        method_names[1].next_to(circle_group[1]).shift(RIGHT)
        method_names[2].next_to(circle_group[2], direction = LEFT).shift(LEFT)
        method_names[3].next_to(circle_group[3], direction = LEFT).shift(LEFT)
        method_names[4].next_to(circle_group[4]).shift(RIGHT)

        self.add(circle_group, method_names)
        center_dot = Dot(circle_group[3].get_center(), color = GREEN, stroke_color = BLUE, stroke_width = 3)
        self.play(
            method_names[3].animate.shift(3 * UP + 0.25 * RIGHT),
            Write(center_dot)
        )
        self.wait()

        def get_chord_by_randompointrandomangle(circle, randomness = True):
            chord_group = VGroup()
            circle_center = circle.get_center()
            circle_radius = circle.radius
            r = circle_radius * np.sqrt(np.random.random()) if randomness else rad / 2
            semichord = np.sqrt(circle_radius ** 2 - r * r)
            first_point = circle_center + r * RIGHT + semichord * UP
            second_point = circle_center + r * RIGHT + semichord * DOWN
            perpendicular_chord = Line(first_point, second_point, stroke_width = 0.875)
            chord_group.add(perpendicular_chord)
            horizontal_line = Line(circle_group[3].get_center(), circle_group[3].get_center() + rad * RIGHT, color = GREEN, stroke_width = 0.875)
            chord_group.add(horizontal_line)
            random_angle = 2 * PI * np.random.random() if randomness else 45 * DEGREES
            r = r * np.cos(random_angle)
            semichord = np.sqrt(circle_radius ** 2 - r * r)
            first_point = circle_center + r * RIGHT + semichord * UP
            second_point = circle_center + r * RIGHT + semichord * DOWN
            chord = Line(first_point, second_point).rotate(random_angle, about_point = circle_group[3].get_center())
            if chord.get_length() >= rad * np.sqrt(3):
                chord.color = BLUE
            else:
                chord.color = WHITE
            chord_group.add(chord)
            chord_group.rotate(2 * PI * np.random.random(), about_point = circle_group[3].get_center())
            return chord_group
        
        n = 1000
        chord_groups = VGroup(*[get_chord_by_randompointrandomangle(circle_group[3]) for _ in range(n)])
        chords, perpendiculars, horizontals = VGroup(), VGroup(), VGroup()
        for obj in chord_groups:
            perpendiculars.add(obj[0])
            horizontals.add(obj[1])
            chords.add(obj[-1])
        
        illustrative_chord = get_chord_by_randompointrandomangle(circle_group[3], randomness = False)
        glow = VGroup(*[Line(ORIGIN, 1 * RIGHT / 4, color = YELLOW, stroke_width = 1 / 4).rotate_about_origin(k * 360 * DEGREES / 100) for k in range(100)])
        illustrative_glow = glow.copy().move_to(illustrative_chord[0].get_center())

        self.play(FadeIn(illustrative_glow))
        self.wait()
        self.play(Write(illustrative_chord[:2]))
        self.play(Write(illustrative_chord[2]))
        self.wait()

        thicker_chords = VGroup()
        for cho, perpcho, radl in zip(chords, perpendiculars, horizontals):
            temp_group = VGroup()
            temp_group.add(cho.copy())
            cho.stroke_width = 0.3125
            cho.stroke_opacity = 0.25
            first_glow = glow.copy().move_to(perpcho.get_center())
            temp_group.add(first_glow, perpcho, radl)
            thicker_chords.add(temp_group)
        
        self.play(FadeOut(VGroup(illustrative_chord, illustrative_glow)), run_time = 0.1)
        self.play(
            ShowIncreasingSubsets(chords),
            ShowSubmobjectsOneByOne(thicker_chords),
            run_time = 10,
            rate_func = linear
        )
        self.play(FadeOut(thicker_chords), run_time = 0.1)
        self.wait()
        illustrative_angle = illustrative_chord[1].get_angle()
        illustrative_chord.rotate(-illustrative_angle, about_point = circle_group[3].get_center())
        illustrative_glow.rotate(-illustrative_angle, about_point = circle_group[3].get_center())
        self.play(FadeIn(illustrative_chord), FadeIn(illustrative_glow), FadeOut(chords))
        self.wait()

        x = rad / 2
        y = np.sqrt(rad ** 2 - x ** 2)
        horizontal_brace = Brace(Line(ORIGIN, x * RIGHT)).shift(circle_group[3].get_center())
        smallx = MathTex("r")
        smallx.set_background_stroke(width = 5)
        smallx.next_to(horizontal_brace, DOWN)
        chord_length = MathTex("2\\sqrt{1-(r\\cos t)^2}")
        chord_length.set_background_stroke(width = 5)
        chord_length.rotate(- PI / 4)
        chord_length.move_to(illustrative_chord[-1])
        chord_length.shift(1 * UR / 2 / np.sqrt(2))
        illustrative_t = MathTex("t")
        illustrative_t.set_background_stroke(width = 5)
        illustrative_tangle = Angle(illustrative_chord[0], illustrative_chord[-1], radius = 0.5)
        illustrative_t.next_to(illustrative_tangle, direction = DOWN)
        expected_length = VGroup(
            MathTex("\\mathbb{E}(L)", "=", "\\int_0^{\\pi/2}\\int_0^1 2\\sqrt{1-(r\\cos t)^2}\\cdot 2r\\,dr\\frac{dt}{\\pi/2}"),
            MathTex("=", "\\frac{16}{3\\pi} \\approx 1.6977")
        )
        expected_length.move_to(circle_group[3].get_center() + 6 * UP / 4 + 7 * LEFT)
        for obj in expected_length:
            obj.scale(0.875)
        expected_length[0].shift(RIGHT + RIGHT / 4)
        expected_length[1].move_to(expected_length[0][1].get_center() + expected_length[1].get_center() - expected_length[1][0].get_center())
        expected_length[1].shift(1.5 * DOWN)
        expected_length[0].add_background_rectangle()

        self.play(
            FadeIn(VGroup(horizontal_brace, illustrative_tangle)),
            FadeIn(smallx), Write(chord_length), Write(illustrative_t),
            Write(expected_length[0]),
            FadeIn(chords)
        )
        self.wait()
        self.play(Write(expected_length[1]))
        self.wait()

        vertical_chord = illustrative_chord[-1].copy()
        self.play(
            FadeOut(VGroup(illustrative_glow, horizontal_brace, illustrative_tangle, chord_length, smallx, illustrative_t, illustrative_chord)),
            FadeIn(vertical_chord)
        )
        self.wait()

        animlist = []
        for obj in chords:
            angle = Line(circle_group[3].get_center(), obj.get_end()).get_angle()
            animlist += [Rotate(obj, -angle - PI / 2, about_point = circle_group[3].get_center())]
        angle = Line(circle_group[3].get_center(), vertical_chord.get_end()).get_angle()
        self.play(Rotate(vertical_chord, -angle - PI / 2, about_point = circle_group[3].get_center()))
        self.play(LaggedStart(*animlist), rate_func = linear, run_time = 5)
        self.wait()

        tangent_line = Line(vertical_chord.get_end() + 2 * LEFT, vertical_chord.get_end() + 2 * RIGHT, color = GREEN)
        center_line = Line(circle_group[3].get_center() + rad * DOWN + DOWN / 4, circle_group[3].get_center() + rad * UP + UP / 4, color = RED)
        tangent_angle = Angle(tangent_line, Line(vertical_chord.get_end(), vertical_chord.get_start()), radius = 0.75)
        theta_text = MathTex("\\theta").set_background_stroke(width = 5)
        theta_text.next_to(tangent_angle)
        theta_text.shift(UP / 4)
        theta_distribution = MathTex("f_{\\Theta}(\\theta)", "=", "\\frac{2\\sin^2\\theta}{\\pi/2}")
        theta_distribution.move_to(expected_length[-1])
        theta_distribution.shift(2 * DOWN)
        self.play(
            Write(theta_distribution),
            Write(theta_text),
            FadeIn(tangent_line, center_line, tangent_angle)
        )
        self.wait()

        self.play(
            FadeOut(VGroup(tangent_line, center_line, center_dot, tangent_angle, theta_text, chords, vertical_chord)),
            self.camera.frame.animate.move_to(circle_group[4].get_center() + 3 * RIGHT)
        )
        self.wait()

class RandomRadialRandomAngle2(MovingCameraScene):
    def construct(self):
        rad = 3.75
        pentagon_rad = 7
        base_circle = Circle(radius = rad, color = WHITE).shift(pentagon_rad * RIGHT)
        circle_group = VGroup()
        for k in range(5):
            circle_group.add(base_circle.copy().rotate_about_origin(k * 72 * DEGREES))
        self.camera.frame.move_to(circle_group[4].get_center() + 3 * RIGHT)

        method_names = VGroup(
            MathTex("\\text{Circumference method}"),
            MathTex("\\text{Radial line method}"),
            MathTex("\\text{Midpoint method}"),
            MathTex("\\text{Random point on circle,}\\\\", "\\text{line at random angle}"),
            MathTex("\\text{Random point on random radii,}\\\\", "\\text{line at random angle}"),
        )

        method_names[0].next_to(circle_group[0]).shift(RIGHT)
        method_names[1].next_to(circle_group[1]).shift(RIGHT)
        method_names[2].next_to(circle_group[2], direction = LEFT).shift(LEFT)
        method_names[3].next_to(circle_group[3], direction = LEFT).shift(LEFT)
        method_names[4].next_to(circle_group[4]).shift(RIGHT)
        method_names[4].add_background_rectangle()

        self.add(circle_group, method_names)
        center_dot = Dot(circle_group[4].get_center(), color = GREEN, stroke_color = BLUE, stroke_width = 3)
        self.play(
            method_names[4].animate.shift(3 * UP + 1.75 * LEFT),
            Write(center_dot)
        )
        self.wait()

        def get_chord_by_randompointrandomangle(circle, randomness = True):
            chord_group = VGroup()
            circle_center = circle.get_center()
            circle_radius = circle.radius
            r = circle_radius * np.random.random() if randomness else rad / 2
            semichord = np.sqrt(circle_radius ** 2 - r * r)
            first_point = circle_center + r * RIGHT + semichord * UP
            second_point = circle_center + r * RIGHT + semichord * DOWN
            perpendicular_chord = Line(first_point, second_point, stroke_width = 0.875)
            chord_group.add(perpendicular_chord)
            horizontal_line = Line(circle_group[4].get_center(), circle_group[4].get_center() + rad * RIGHT, color = GREEN, stroke_width = 0.875)
            chord_group.add(horizontal_line)
            random_angle = 2 * PI * np.random.random() if randomness else 45 * DEGREES
            r = r * np.cos(random_angle)
            semichord = np.sqrt(circle_radius ** 2 - r * r)
            first_point = circle_center + r * RIGHT + semichord * UP
            second_point = circle_center + r * RIGHT + semichord * DOWN
            chord = Line(first_point, second_point).rotate(random_angle, about_point = circle_group[4].get_center())
            if chord.get_length() >= rad * np.sqrt(3):
                chord.color = BLUE
            else:
                chord.color = WHITE
            chord_group.add(chord)
            chord_group.rotate(2 * PI * np.random.random(), about_point = circle_group[4].get_center())
            return chord_group
        
        n = 1000
        chord_groups = VGroup(*[get_chord_by_randompointrandomangle(circle_group[4]) for _ in range(n)])
        chords, perpendiculars, horizontals = VGroup(), VGroup(), VGroup()
        for obj in chord_groups:
            perpendiculars.add(obj[0])
            horizontals.add(obj[1])
            chords.add(obj[-1])
        
        illustrative_chord = get_chord_by_randompointrandomangle(circle_group[4], randomness = False)
        glow = VGroup(*[Line(ORIGIN, 1 * RIGHT / 4, color = YELLOW, stroke_width = 1 / 4).rotate_about_origin(k * 360 * DEGREES / 100) for k in range(100)])
        illustrative_glow = glow.copy().move_to(illustrative_chord[0].get_center())

        self.play(FadeIn(illustrative_glow))
        self.wait()
        self.play(Write(illustrative_chord[:2]))
        self.play(Write(illustrative_chord[2]))
        self.wait()

        thicker_chords = VGroup()
        for cho, perpcho, radl in zip(chords, perpendiculars, horizontals):
            temp_group = VGroup()
            temp_group.add(cho.copy())
            cho.stroke_width = 0.3125
            cho.stroke_opacity = 0.25
            first_glow = glow.copy().move_to(perpcho.get_center())
            temp_group.add(first_glow, perpcho, radl)
            thicker_chords.add(temp_group)
        
        self.play(FadeOut(VGroup(illustrative_chord, illustrative_glow)), run_time = 0.1)
        self.play(
            ShowIncreasingSubsets(chords),
            ShowSubmobjectsOneByOne(thicker_chords),
            run_time = 10,
            rate_func = linear
        )
        self.play(FadeOut(thicker_chords), run_time = 0.1)
        self.wait()
        illustrative_angle = illustrative_chord[1].get_angle()
        illustrative_chord.rotate(-illustrative_angle + PI, about_point = circle_group[4].get_center())
        illustrative_glow.rotate(-illustrative_angle + PI, about_point = circle_group[4].get_center())
        self.play(FadeIn(illustrative_chord), FadeIn(illustrative_glow), FadeOut(chords))
        self.wait()

        x = rad / 2
        y = np.sqrt(rad ** 2 - x ** 2)
        horizontal_brace = Brace(Line(ORIGIN, x * LEFT), direction = UP).shift(circle_group[4].get_center())
        smallx = MathTex("x")
        smallx.set_background_stroke(width = 5)
        smallx.next_to(horizontal_brace, UP)
        chord_length = MathTex("2\\sqrt{1-(x\\cos t)^2}")
        chord_length.set_background_stroke(width = 5)
        chord_length.rotate(- PI / 4)
        chord_length.move_to(illustrative_chord[-1])
        chord_length.shift(1 * DL / 2 / np.sqrt(2))
        illustrative_t = MathTex("t")
        illustrative_t.set_background_stroke(width = 5)
        illustrative_tangle = Angle(illustrative_chord[0], illustrative_chord[-1], radius = 0.5)
        illustrative_t.next_to(illustrative_tangle, direction = UP)
        expected_length = VGroup(
            MathTex("\\mathbb{E}(L)", "=", "\\int_0^{\\pi/2}\\int_0^1 2\\sqrt{1-(x\\cos t)^2}\\cdot \\,dx\\frac{dt}{\\pi/2}"),
            MathTex("=", "\\frac{2}{\\pi}(2G+1) \\approx 1.8029")
        )
        expected_length.move_to(circle_group[4].get_center() + 6 * UP / 4 + 7 * RIGHT)
        for obj in expected_length:
            obj.scale(0.875)
        expected_length[0].shift(LEFT)
        expected_length[1].move_to(expected_length[0][1].get_center() + expected_length[1].get_center() - expected_length[1][0].get_center())
        expected_length[1].shift(1.5 * DOWN + 1.5 * RIGHT)
        expected_length[0].add_background_rectangle()

        self.play(
            FadeIn(VGroup(horizontal_brace, illustrative_tangle)),
            Write(VGroup(chord_length, smallx, illustrative_t)),
            Write(expected_length[0]),
            FadeIn(chords)
        )
        self.wait()
        self.play(Write(expected_length[1]))
        self.wait()

        vertical_chord = illustrative_chord[-1].copy()
        self.play(
            FadeOut(VGroup(illustrative_glow, horizontal_brace, illustrative_tangle, chord_length, smallx, illustrative_t, illustrative_chord)),
            FadeIn(vertical_chord)
        )
        self.wait()

        animlist = []
        for obj in chords:
            angle = Line(circle_group[4].get_center(), obj.get_end()).get_angle()
            animlist += [Rotate(obj, -angle - PI / 2, about_point = circle_group[4].get_center())]
        angle = Line(circle_group[4].get_center(), vertical_chord.get_end()).get_angle()
        self.play(Rotate(vertical_chord, -angle - PI / 2, about_point = circle_group[4].get_center()))
        self.play(LaggedStart(*animlist), rate_func = linear, run_time = 5)
        self.wait()

        tangent_line = Line(vertical_chord.get_end() + 2 * LEFT, vertical_chord.get_end() + 2 * RIGHT, color = GREEN)
        center_line = Line(circle_group[4].get_center() + rad * DOWN + DOWN / 4, circle_group[4].get_center() + rad * UP + UP / 4, color = RED)
        tangent_angle = Angle(tangent_line, Line(vertical_chord.get_end(), vertical_chord.get_start()), radius = 0.75)
        theta_text = MathTex("\\theta").set_background_stroke(width = 5)
        theta_text.next_to(tangent_angle)
        theta_text.shift(UP / 4)
        theta_distribution = MathTex("f_{\\Theta}(\\theta)", "=", "\\frac{\\sin\\theta\\sinh^{-1}(\\tan\\theta)}{\\pi/2}")
        theta_distribution.move_to(expected_length[-1])
        theta_distribution.shift(2 * DOWN)
        self.play(
            Write(theta_distribution),
            Write(theta_text),
            FadeIn(tangent_line, center_line, tangent_angle)
        )
        self.wait()

        self.play(
            FadeOut(VGroup(tangent_line, center_line, center_dot, tangent_angle, theta_text, chords, vertical_chord)),
            self.camera.frame.animate.move_to(circle_group[0].get_center() + 3 * RIGHT),
            FadeOut(method_names[4])
        )
        self.wait()

class CircumferenceMethodProbability1(MovingCameraScene):
    def construct(self):
        rad = 3.75
        pentagon_rad = 7
        base_circle = Circle(radius = rad, color = WHITE).shift(pentagon_rad * RIGHT)
        circle_group = VGroup()
        for k in range(5):
            circle_group.add(base_circle.copy().rotate_about_origin(k * 72 * DEGREES))
        self.camera.frame.scale(2)
        self.play(Write(circle_group))
        self.wait()

        method_names = VGroup(
            MathTex("\\text{Circumference method}"),
            MathTex("\\text{Radial line method}"),
            MathTex("\\text{Midpoint method}"),
            MathTex("\\text{Random point on circle,}\\\\", "\\text{line at random angle}"),
            MathTex("\\text{Random point on random radii,}\\\\", "\\text{line at random angle}"),
        )

        method_names[0].next_to(circle_group[0]).shift(RIGHT)
        method_names[1].next_to(circle_group[1]).shift(RIGHT)
        method_names[2].next_to(circle_group[2], direction = LEFT).shift(LEFT)
        method_names[3].next_to(circle_group[3], direction = LEFT).shift(LEFT)
        method_names[4].next_to(circle_group[4]).shift(RIGHT)

        self.play(
            self.camera.frame.animate.move_to(circle_group[0].get_center() + 3 * RIGHT).scale(1 / 2),
            #Write(method_names)
        )
        method_names[0].to_edge(UP).shift(0.5 * LEFT)
        temp_dot = Dot(circle_group[0].get_right())
        temp_dot.rotate(45 * DEGREES, about_point = circle_group[0].get_center())
        random_chord = Line(circle_group[0].get_bottom(), temp_dot.get_center(), color = BLUE)
        def add_end_dots(line):
            obj = line.copy()
            startingdot = Dot(obj.get_start(), color = YELLOW, stroke_color = RED, stroke_width = 3)
            endingdot = Dot(obj.get_end(), color = YELLOW, stroke_color = RED, stroke_width = 3)
            return VGroup(obj, startingdot, endingdot)
        random_chord = add_end_dots(random_chord)
        self.play(Write(random_chord))
        self.wait()

        event_definition = VGroup(
            MathTex("E:\\text{points chosen }", "\\text{uniformly randomly}"),
            MathTex("\\text{on a circle lies on different sides of a}"),
            MathTex("\\text{random}", "\\text{ chord}")
        )
        event_definition.arrange(DOWN)
        event_definition[-1].shift(DOWN / 8)
        event_definition.move_to(circle_group[0].get_center() + 10 * UP / 4 + 5.5 * RIGHT)
        event_definition.add_background_rectangle()
        self.play(Write(event_definition))
        self.play(
            event_definition[1][-1].animate.set_color(GREEN),
            event_definition[3][0].animate.set_color(YELLOW)
        )
        self.wait()

        tangent_line = Line(circle_group[0].get_bottom() + LEFT, circle_group[0].get_bottom() + RIGHT, color = TEAL)
        tangent_angle = Angle(tangent_line, random_chord[0], radius = 0.5)
        theta_text = MathTex("\\theta").set_background_stroke(width = 5)
        theta_text.next_to(tangent_angle)
        theta_text.shift(UP / 4)
        area_points = [
            [0.6654, -3.1073, 0],
            [1.712, -3.0153, 0],
            [1, -2, 0],
            [2.5748, -2.0146, 0],
            [1.6316, -0.7723, 0],
            [3.3569, -0.6918, 0],
            [2.0917, 0.4009, 0],
            [3.3339, 0.6539, 0],
            [2.5057, 1.4591, 0],
            [3, 2, 0]
        ]
        area_marker = VMobject(color = RED, stroke_width = 10).set_points_smoothly(area_points)
        area_marker.shift(circle_group[0].get_center())
        self.play(Create(area_marker), run_time = 2, rate_func = linear)
        self.play(
            Write(VGroup(tangent_angle, tangent_line, theta_text))
        )
        self.wait()

        probability_text = MathTex("\\text{Area}(\\theta)", "=", "\\frac{2\\theta-\\sin 2\\theta}{2}")
        event_probability = MathTex("\\mathbb{P}(E|\\theta)", "=", "2\\left(\\frac{2\\theta-\\sin 2 \\theta}{2\\pi}\\right)\\left(1-\\frac{2\\theta-\\sin2\\theta}{2\\pi}\\right)")
        probability_text.move_to(circle_group[0].get_center() + 1 * UP / 4 + 7 * RIGHT)
        event_probability.move_to(circle_group[0].get_center() - 8 * UP / 4 + 4.5 * RIGHT)
        event_probability.add_background_rectangle()

        self.play(Write(probability_text))
        self.wait()
        self.play(Write(event_probability))
        self.wait()

        self.play(
            FadeOut(area_marker),
            FadeOut(VGroup(event_probability, probability_text)),
            #event_probability.animate.move_to(circle_group[0].get_center() + 12 * UP / 4 + 5 * RIGHT)
        )
        self.wait()

        green_glow = VGroup(*[Line(ORIGIN, 1 * RIGHT / 4, color = GREEN, stroke_width = 1 / 4).rotate_about_origin(k * 360 * DEGREES / 100) for k in range(100)])
        red_glow = VGroup(*[Line(ORIGIN, 1 * RIGHT / 4, color = RED, stroke_width = 1 / 4).rotate_about_origin(k * 360 * DEGREES / 100) for k in range(100)])
        def get_two_points(circle):
            circle_center = circle.get_center()
            circle_radius = circle.radius
            temp_group = VGroup()
            temp_dot1 = Dot(circle_center + RIGHT * circle_radius * np.sqrt(np.random.random())).rotate(2 * PI * np.random.random(), about_point = circle_center)
            temp_dot2 = Dot(circle_center + RIGHT * circle_radius * np.sqrt(np.random.random())).rotate(2 * PI * np.random.random(), about_point = circle_center)
            temp_dot1_lies_right =  temp_dot1.get_center()[0] >= circle_center[0] + 1.4350628
            temp_dot2_lies_right =  temp_dot2.get_center()[0] >= circle_center[0] + 1.4350628
            good_dots = not (temp_dot1_lies_right ^ temp_dot2_lies_right)
            temp_group.add(temp_dot1, temp_dot2)
            if good_dots:
                temp_group.add(red_glow.copy().move_to(temp_dot1), red_glow.copy().move_to(temp_dot2))
            else:
                temp_group.add(green_glow.copy().move_to(temp_dot1), green_glow.copy().move_to(temp_dot2))
            temp_group.rotate(-22.5 * DEGREES, about_point = circle_center)
            return temp_group
        n = 50
        points_group = VGroup(*[get_two_points(circle_group[0]) for _ in range(n)])
        self.play(
            ShowSubmobjectsOneByOne(points_group),
            run_time = 10,
            rate_func = linear
        )
        self.wait()

        density_function = MathTex("f(\\theta)", "=", "\\frac{1}{\\pi/2}")
        probability_statements = VGroup(
            MathTex("\\mathbb{P}(E)", "=", "\\int_0^{\\pi/2}\\mathbb{P}(E|\\theta)\\,\\frac{d\\theta}{\\pi/2}"),
            MathTex("=", "\\frac{1}{3}-\\frac{5}{4\\pi^2}"),
            #MathTex("\\approx", "0.2067")
        )
        density_function.move_to(circle_group[0].get_center() + 2 * UP / 4 + 6.75 * RIGHT)
        box = SurroundingRectangle(density_function)
        probability_statements[0].move_to(circle_group[0].get_center() - 5 * UP / 4 + 6.75 * RIGHT)
        probability_statements[1].move_to(probability_statements[0][1].get_center() + probability_statements[1].get_center() - probability_statements[1][0].get_center())
        #probability_statements[2].move_to(probability_statements[0][1].get_center() + probability_statements[2].get_center() - probability_statements[2][0].get_center())
        probability_statements[1].shift(3 * DOWN / 2)
        #probability_statements[2].shift(6 * DOWN / 2)

        self.play(Write(density_function), Create(box))
        self.wait()
        self.play(Write(probability_statements[0]))
        self.wait()
        self.play(Write(probability_statements[1:]))
        self.wait()
        self.play(
            self.camera.frame.animate.move_to(circle_group[1].get_center() + 3 * RIGHT),
            event_definition.animate.move_to(circle_group[1].get_center() + 10 * UP / 4 + 5.5 * RIGHT),
            FadeOut(VGroup(random_chord, points_group))
        )
        self.wait()

class RadialLineProbability(MovingCameraScene):
    def construct(self):
        rad = 3.75
        pentagon_rad = 7
        base_circle = Circle(radius = rad, color = WHITE).shift(pentagon_rad * RIGHT)
        circle_group = VGroup()
        for k in range(5):
            circle_group.add(base_circle.copy().rotate_about_origin(k * 72 * DEGREES))
        self.camera.frame.scale(2)
        self.add(circle_group)

        method_names = VGroup(
            MathTex("\\text{Circumference method}"),
            MathTex("\\text{Radial line method}"),
            MathTex("\\text{Midpoint method}"),
            MathTex("\\text{Random point on circle,}\\\\", "\\text{line at random angle}"),
            MathTex("\\text{Random point on random radii,}\\\\", "\\text{line at random angle}"),
        )

        self.camera.frame.move_to(circle_group[1].get_center() + 3 * RIGHT).scale(1 / 2)
        method_names[1].move_to(circle_group[1].get_center() + 3.75 * DOWN + 2 * LEFT)
        #method_names[1].move_to(circle_group[1].get_center())
        method_names[1].add_background_rectangle()
        temp_dot = Dot(circle_group[1].get_center() + rad * RIGHT)
        temp_dot.rotate(45 * DEGREES, about_point = circle_group[1].get_center())
        random_chord = Line(circle_group[1].get_center() + rad * DOWN, temp_dot.get_center(), color = BLUE)
        event_definition = VGroup(
            MathTex("E:\\text{points chosen }", "\\text{uniformly randomly}"),
            MathTex("\\text{on a circle lies on different sides of a}"),
            MathTex("\\text{random}", "\\text{ chord}")
        )
        event_definition.arrange(DOWN)
        event_definition[-1].shift(DOWN / 8)
        event_definition.move_to(circle_group[1].get_center() + 10 * UP / 4 + 5.5 * RIGHT)
        event_definition.add_background_rectangle()
        event_definition[1][-1].set_color(GREEN),
        event_definition[3][0].set_color(YELLOW)

        def add_end_dots(line):
            obj = line.copy()
            startingdot = Dot(obj.get_start(), color = YELLOW, stroke_color = RED, stroke_width = 3)
            endingdot = Dot(obj.get_end(), color = YELLOW, stroke_color = RED, stroke_width = 3)
            return VGroup(obj, startingdot, endingdot)

        random_chord = add_end_dots(random_chord)
        tangent_line = Line(circle_group[1].get_center() + rad * DOWN + LEFT, circle_group[1].get_center() + rad * DOWN + RIGHT, color = TEAL)
        tangent_angle = Angle(tangent_line, random_chord[0], radius = 0.5)
        theta_text = MathTex("\\theta").set_background_stroke(width = 5)
        theta_text.next_to(tangent_angle)
        theta_text.shift(UP / 4)
        
        self.add(event_definition)
        self.add_foreground_mobject(event_definition)
        self.wait()
        self.play(
            Write(random_chord), Write(theta_text),
            FadeIn(VGroup(tangent_angle, tangent_line)),
            FadeIn(method_names[1])
        )
        self.wait()

        green_glow = VGroup(*[Line(ORIGIN, 1 * RIGHT / 4, color = GREEN, stroke_width = 1 / 4).rotate_about_origin(k * 360 * DEGREES / 100) for k in range(100)])
        red_glow = VGroup(*[Line(ORIGIN, 1 * RIGHT / 4, color = RED, stroke_width = 1 / 4).rotate_about_origin(k * 360 * DEGREES / 100) for k in range(100)])
        def get_two_points(circle):
            circle_center = circle.get_center()
            circle_radius = circle.radius
            temp_group = VGroup()
            temp_dot1 = Dot(circle_center + RIGHT * circle_radius * np.sqrt(np.random.random())).rotate(2 * PI * np.random.random(), about_point = circle_center)
            temp_dot2 = Dot(circle_center + RIGHT * circle_radius * np.sqrt(np.random.random())).rotate(2 * PI * np.random.random(), about_point = circle_center)
            temp_dot1_lies_right =  temp_dot1.get_center()[0] >= circle_center[0] + 1.4350628
            temp_dot2_lies_right =  temp_dot2.get_center()[0] >= circle_center[0] + 1.4350628
            good_dots = not (temp_dot1_lies_right ^ temp_dot2_lies_right)
            temp_group.add(temp_dot1, temp_dot2)
            if good_dots:
                temp_group.add(red_glow.copy().move_to(temp_dot1), red_glow.copy().move_to(temp_dot2))
            else:
                temp_group.add(green_glow.copy().move_to(temp_dot1), green_glow.copy().move_to(temp_dot2))
            temp_group.rotate(-22.5 * DEGREES, about_point = circle_center)
            return temp_group
        n = 50
        points_group = VGroup(*[get_two_points(circle_group[1]) for _ in range(n)])
        self.play(
            ShowSubmobjectsOneByOne(points_group),
            run_time = 10,
            rate_func = linear
        )
        self.wait()

        density_function = MathTex("f(\\theta)", "=", "\\sin\\theta")
        probability_statements = VGroup(
            MathTex("\\mathbb{P}(E)", "=", "\\int_0^{\\pi/2}\\mathbb{P}(E|\\theta)f(\\theta)\\,d\\theta"),
            MathTex("=", "\\frac{128}{45\\pi^2}"),
            MathTex("\\approx", "0.2882")
        )
        density_function.move_to(circle_group[1].get_center() + 9 * UP / 4 + 6.75 * RIGHT)
        box = SurroundingRectangle(density_function)
        #probability_statements[0].set_background_stroke(width = 5)
        probability_statements[0].move_to(circle_group[1].get_center() + 2 * UP / 4 + 6.75 * RIGHT)
        probability_statements[1].move_to(probability_statements[0][1].get_center() + probability_statements[1].get_center() - probability_statements[1][0].get_center())
        probability_statements[2].move_to(probability_statements[0][1].get_center() + probability_statements[2].get_center() - probability_statements[2][0].get_center())
        probability_statements[1].shift(3 * DOWN / 2)
        probability_statements[2].shift(6 * DOWN / 2)

        self.play(
            Write(density_function),
            Create(box),
            event_definition.animate.shift(5 * LEFT)
        )
        self.wait()
        self.play(Write(probability_statements[0]))
        self.wait()
        self.play(Write(probability_statements[1:]))
        self.wait()
        self.play(
            self.camera.frame.animate.move_to(circle_group[2].get_center() + 3 * LEFT),
            event_definition.animate.move_to(circle_group[2].get_center() + 10 * UP / 4 - 5.5 * RIGHT),
            FadeOut(VGroup(random_chord, points_group, theta_text, tangent_angle, tangent_line, method_names[1]))
        )
        self.wait()

class MidpointProbability(MovingCameraScene):
    def construct(self):
        rad = 3.75
        pentagon_rad = 7
        base_circle = Circle(radius = rad, color = WHITE).shift(pentagon_rad * RIGHT)
        circle_group = VGroup()
        for k in range(5):
            circle_group.add(base_circle.copy().rotate_about_origin(k * 72 * DEGREES))
        self.camera.frame.scale(2)
        self.add(circle_group)

        method_names = VGroup(
            MathTex("\\text{Circumference method}"),
            MathTex("\\text{Radial line method}"),
            MathTex("\\text{Midpoint method}"),
            MathTex("\\text{Random point on circle,}\\\\", "\\text{line at random angle}"),
            MathTex("\\text{Random point on random radii,}\\\\", "\\text{line at random angle}"),
        )

        j = 2

        self.camera.frame.move_to(circle_group[j].get_center() + 3 * LEFT).scale(1 / 2)
        method_names[j].move_to(circle_group[j].get_center() + 3.75 * DOWN + 2 * RIGHT)
        method_names[j].set_background_stroke(width = 5)
        temp_dot = Dot(circle_group[j].get_center() + rad * RIGHT)
        temp_dot.rotate(45 * DEGREES, about_point = circle_group[j].get_center())
        random_chord = Line(circle_group[j].get_center() + rad * DOWN, temp_dot.get_center(), color = BLUE)
        event_definition = VGroup(
            MathTex("E:\\text{points chosen }", "\\text{uniformly randomly}"),
            MathTex("\\text{on a circle lies on different sides of a}"),
            MathTex("\\text{random}", "\\text{ chord}")
        )
        event_definition.arrange(DOWN)
        event_definition[-1].shift(DOWN / 8)
        event_definition.move_to(circle_group[j].get_center() + 10 * UP / 4 + 5.5 * LEFT)
        event_definition.add_background_rectangle()
        event_definition[1][-1].set_color(GREEN),
        event_definition[3][0].set_color(YELLOW)

        def add_end_dots(line):
            obj = line.copy()
            startingdot = Dot(obj.get_start(), color = YELLOW, stroke_color = RED, stroke_width = 3)
            endingdot = Dot(obj.get_end(), color = YELLOW, stroke_color = RED, stroke_width = 3)
            return VGroup(obj, startingdot, endingdot)

        random_chord = add_end_dots(random_chord)
        tangent_line = Line(circle_group[j].get_center() + rad * DOWN + LEFT, circle_group[j].get_center() + rad * DOWN + RIGHT, color = TEAL)
        tangent_angle = Angle(tangent_line, random_chord[0], radius = 0.5)
        theta_text = MathTex("\\theta").set_background_stroke(width = 5)
        theta_text.next_to(tangent_angle)
        theta_text.shift(UP / 4)
        
        self.add(event_definition)
        self.add_foreground_mobject(event_definition)
        self.wait()
        self.play(
            Write(random_chord), Write(theta_text),
            FadeIn(VGroup(tangent_angle, tangent_line)),
            FadeIn(method_names[j])
        )
        self.wait()

        green_glow = VGroup(*[Line(ORIGIN, 1 * RIGHT / 4, color = GREEN, stroke_width = 1 / 4).rotate_about_origin(k * 360 * DEGREES / 100) for k in range(100)])
        red_glow = VGroup(*[Line(ORIGIN, 1 * RIGHT / 4, color = RED, stroke_width = 1 / 4).rotate_about_origin(k * 360 * DEGREES / 100) for k in range(100)])
        def get_two_points(circle):
            circle_center = circle.get_center()
            circle_radius = circle.radius
            temp_group = VGroup()
            temp_dot1 = Dot(circle_center + RIGHT * circle_radius * np.sqrt(np.random.random())).rotate(2 * PI * np.random.random(), about_point = circle_center)
            temp_dot2 = Dot(circle_center + RIGHT * circle_radius * np.sqrt(np.random.random())).rotate(2 * PI * np.random.random(), about_point = circle_center)
            temp_dot1_lies_right =  temp_dot1.get_center()[0] >= circle_center[0] + 1.4350628
            temp_dot2_lies_right =  temp_dot2.get_center()[0] >= circle_center[0] + 1.4350628
            good_dots = not (temp_dot1_lies_right ^ temp_dot2_lies_right)
            temp_group.add(temp_dot1, temp_dot2)
            if good_dots:
                temp_group.add(red_glow.copy().move_to(temp_dot1), red_glow.copy().move_to(temp_dot2))
            else:
                temp_group.add(green_glow.copy().move_to(temp_dot1), green_glow.copy().move_to(temp_dot2))
            temp_group.rotate(-22.5 * DEGREES, about_point = circle_center)
            return temp_group
        n = 50
        points_group = VGroup(*[get_two_points(circle_group[j]) for _ in range(n)])
        self.play(
            ShowSubmobjectsOneByOne(points_group),
            run_time = 10,
            rate_func = linear
        )
        self.wait()

        density_function = MathTex("f(\\theta)", "=", "\\sin 2\\theta")
        probability_statements = VGroup(
            MathTex("\\mathbb{P}(E)", "=", "\\int_0^{\\pi/2}\\mathbb{P}(E|\\theta)f(\\theta)\\,d\\theta"),
            MathTex("=", "\\frac{1}{8} + \\frac{2}{3\\pi^2}"),
            MathTex("\\approx", "0.1926")
        )
        density_function.move_to(circle_group[j].get_center() + 9 * UP / 4 + 6.75 * LEFT)
        box = SurroundingRectangle(density_function)
        probability_statements[0].set_background_stroke(width = 5)
        probability_statements[0].move_to(circle_group[j].get_center() + 2 * UP / 4 + 6.75 * LEFT)
        probability_statements[1].move_to(probability_statements[0][1].get_center() + probability_statements[1].get_center() - probability_statements[1][0].get_center())
        probability_statements[2].move_to(probability_statements[0][1].get_center() + probability_statements[2].get_center() - probability_statements[2][0].get_center())
        probability_statements[1].shift(3 * DOWN / 2)
        probability_statements[2].shift(6 * DOWN / 2)

        self.play(
            Write(density_function),
            Create(box),
            event_definition.animate.shift(5 * RIGHT)
        )
        self.wait()
        self.play(Write(probability_statements[0]))
        self.wait()
        self.play(Write(probability_statements[1:]))
        self.wait()
        self.play(
            self.camera.frame.animate.move_to(circle_group[j + 1].get_center() + 3 * LEFT),
            #event_definition.animate.move_to(circle_group[j + 1].get_center() + 11 * UP / 4 - 5.5 * RIGHT),
            FadeOut(VGroup(event_definition, random_chord, points_group, theta_text, tangent_angle, tangent_line, method_names[1]))
        )
        self.wait()

class RadomPointRandomAngleProbability(MovingCameraScene):
    def construct(self):
        rad = 3.75
        pentagon_rad = 7
        base_circle = Circle(radius = rad, color = WHITE).shift(pentagon_rad * RIGHT)
        circle_group = VGroup()
        for k in range(5):
            circle_group.add(base_circle.copy().rotate_about_origin(k * 72 * DEGREES))
        self.camera.frame.scale(2)
        self.add(circle_group)

        method_names = VGroup(
            MathTex("\\text{Circumference method}"),
            MathTex("\\text{Radial line method}"),
            MathTex("\\text{Midpoint method}"),
            MathTex("\\text{Random point on circle,}\\\\", "\\text{line at random angle}"),
            MathTex("\\text{Random point on random radii,}\\\\", "\\text{line at random angle}"),
        )

        j = 3

        self.camera.frame.move_to(circle_group[j].get_center() + 3 * LEFT).scale(1 / 2)
        method_names[j].move_to(circle_group[j].get_center() + 3.25 * UP + 6 * LEFT)
        method_names[j].set_background_stroke(width = 5)
        temp_dot = Dot(circle_group[j].get_center() + rad * RIGHT)
        temp_dot.rotate(45 * DEGREES, about_point = circle_group[j].get_center())
        random_chord = Line(circle_group[j].get_center() + rad * DOWN, temp_dot.get_center(), color = BLUE)

        def add_end_dots(line):
            obj = line.copy()
            startingdot = Dot(obj.get_start(), color = YELLOW, stroke_color = RED, stroke_width = 3)
            endingdot = Dot(obj.get_end(), color = YELLOW, stroke_color = RED, stroke_width = 3)
            return VGroup(obj, startingdot, endingdot)

        random_chord = add_end_dots(random_chord)
        tangent_line = Line(circle_group[j].get_center() + rad * DOWN + LEFT, circle_group[j].get_center() + rad * DOWN + RIGHT, color = TEAL)
        tangent_angle = Angle(tangent_line, random_chord[0], radius = 0.5)
        theta_text = MathTex("\\theta").set_background_stroke(width = 5)
        theta_text.next_to(tangent_angle)
        theta_text.shift(UP / 4)

        self.play(
            Write(random_chord), Write(theta_text),
            FadeIn(VGroup(tangent_angle, tangent_line)),
            FadeIn(method_names[j])
        )
        self.wait()

        green_glow = VGroup(*[Line(ORIGIN, 1 * RIGHT / 4, color = GREEN, stroke_width = 1 / 4).rotate_about_origin(k * 360 * DEGREES / 100) for k in range(100)])
        red_glow = VGroup(*[Line(ORIGIN, 1 * RIGHT / 4, color = RED, stroke_width = 1 / 4).rotate_about_origin(k * 360 * DEGREES / 100) for k in range(100)])
        def get_two_points(circle):
            circle_center = circle.get_center()
            circle_radius = circle.radius
            temp_group = VGroup()
            temp_dot1 = Dot(circle_center + RIGHT * circle_radius * np.sqrt(np.random.random())).rotate(2 * PI * np.random.random(), about_point = circle_center)
            temp_dot2 = Dot(circle_center + RIGHT * circle_radius * np.sqrt(np.random.random())).rotate(2 * PI * np.random.random(), about_point = circle_center)
            temp_dot1_lies_right =  temp_dot1.get_center()[0] >= circle_center[0] + 1.4350628
            temp_dot2_lies_right =  temp_dot2.get_center()[0] >= circle_center[0] + 1.4350628
            good_dots = not (temp_dot1_lies_right ^ temp_dot2_lies_right)
            temp_group.add(temp_dot1, temp_dot2)
            if good_dots:
                temp_group.add(red_glow.copy().move_to(temp_dot1), red_glow.copy().move_to(temp_dot2))
            else:
                temp_group.add(green_glow.copy().move_to(temp_dot1), green_glow.copy().move_to(temp_dot2))
            temp_group.rotate(-22.5 * DEGREES, about_point = circle_center)
            return temp_group
        n = 50
        points_group = VGroup(*[get_two_points(circle_group[j]) for _ in range(n)])
        self.play(
            ShowSubmobjectsOneByOne(points_group),
            run_time = 10,
            rate_func = linear
        )
        self.wait()

        density_function = MathTex("f(\\theta)", "=", "\\frac{2\\sin^2\\theta}{\\pi/2}")
        probability_statements = VGroup(
            MathTex("\\mathbb{P}(E)", "=", "\\int_0^{\\pi/2}\\mathbb{P}(E|\\theta)f(\\theta)\\,d\\theta"),
            MathTex("=", "\\frac{1}{3}"),
            MathTex("\\approx", "0.333")
        )
        density_function.move_to(circle_group[j].get_center() + 6 * UP / 4 + 6.75 * LEFT)
        box = SurroundingRectangle(density_function)
        probability_statements[0].set_background_stroke(width = 5)
        probability_statements[0].move_to(circle_group[j].get_center() - 1 * UP / 4 + 6.75 * LEFT)
        probability_statements[1].move_to(probability_statements[0][1].get_center() + probability_statements[1].get_center() - probability_statements[1][0].get_center())
        probability_statements[2].move_to(probability_statements[0][1].get_center() + probability_statements[2].get_center() - probability_statements[2][0].get_center())
        probability_statements[1].shift(3 * DOWN / 2)
        probability_statements[2].shift(6 * DOWN / 2)

        self.play(
            Write(density_function),
            Create(box),
            #event_definition.animate.shift(5 * RIGHT)
        )
        self.wait()
        self.play(Write(probability_statements[0]))
        self.wait()
        self.play(Write(probability_statements[1:]))
        self.wait()
        self.play(
            self.camera.frame.animate.move_to(circle_group[j + 1].get_center() + 3 * RIGHT),
            FadeOut(VGroup(random_chord, points_group, theta_text, tangent_angle, tangent_line, method_names[1]))
        )
        self.wait()

class RadomRadialPointRandomAngleProbability(MovingCameraScene):
    def construct(self):
        rad = 3.75
        pentagon_rad = 7
        base_circle = Circle(radius = rad, color = WHITE).shift(pentagon_rad * RIGHT)
        circle_group = VGroup()
        for k in range(5):
            circle_group.add(base_circle.copy().rotate_about_origin(k * 72 * DEGREES))
        self.camera.frame.scale(2)
        self.add(circle_group)

        method_names = VGroup(
            MathTex("\\text{Circumference method}"),
            MathTex("\\text{Radial line method}"),
            MathTex("\\text{Midpoint method}"),
            MathTex("\\text{Random point on circle,}\\\\", "\\text{line at random angle}"),
            MathTex("\\text{Random point on random radii,}\\\\", "\\text{line at random angle}"),
        )

        j = 4

        self.camera.frame.move_to(circle_group[j].get_center() + 3 * RIGHT).scale(1 / 2)
        method_names[j].move_to(circle_group[j].get_center() + 3.25 * UP + 6 * RIGHT)
        method_names[j].set_background_stroke(width = 5)
        temp_dot = Dot(circle_group[j].get_center() + rad * RIGHT)
        temp_dot.rotate(45 * DEGREES, about_point = circle_group[j].get_center())
        random_chord = Line(circle_group[j].get_center() + rad * DOWN, temp_dot.get_center(), color = BLUE)

        def add_end_dots(line):
            obj = line.copy()
            startingdot = Dot(obj.get_start(), color = YELLOW, stroke_color = RED, stroke_width = 3)
            endingdot = Dot(obj.get_end(), color = YELLOW, stroke_color = RED, stroke_width = 3)
            return VGroup(obj, startingdot, endingdot)

        random_chord = add_end_dots(random_chord)
        tangent_line = Line(circle_group[j].get_center() + rad * DOWN + LEFT, circle_group[j].get_center() + rad * DOWN + RIGHT, color = TEAL)
        tangent_angle = Angle(tangent_line, random_chord[0], radius = 0.5)
        theta_text = MathTex("\\theta").set_background_stroke(width = 5)
        theta_text.next_to(tangent_angle)
        theta_text.shift(UP / 4)

        self.play(
            Write(random_chord), Write(theta_text),
            FadeIn(VGroup(tangent_angle, tangent_line)),
            FadeIn(method_names[j])
        )
        self.wait()

        green_glow = VGroup(*[Line(ORIGIN, 1 * RIGHT / 4, color = GREEN, stroke_width = 1 / 4).rotate_about_origin(k * 360 * DEGREES / 100) for k in range(100)])
        red_glow = VGroup(*[Line(ORIGIN, 1 * RIGHT / 4, color = RED, stroke_width = 1 / 4).rotate_about_origin(k * 360 * DEGREES / 100) for k in range(100)])
        def get_two_points(circle):
            circle_center = circle.get_center()
            circle_radius = circle.radius
            temp_group = VGroup()
            temp_dot1 = Dot(circle_center + RIGHT * circle_radius * np.sqrt(np.random.random())).rotate(2 * PI * np.random.random(), about_point = circle_center)
            temp_dot2 = Dot(circle_center + RIGHT * circle_radius * np.sqrt(np.random.random())).rotate(2 * PI * np.random.random(), about_point = circle_center)
            temp_dot1_lies_right =  temp_dot1.get_center()[0] >= circle_center[0] + 1.4350628
            temp_dot2_lies_right =  temp_dot2.get_center()[0] >= circle_center[0] + 1.4350628
            good_dots = not (temp_dot1_lies_right ^ temp_dot2_lies_right)
            temp_group.add(temp_dot1, temp_dot2)
            if good_dots:
                temp_group.add(red_glow.copy().move_to(temp_dot1), red_glow.copy().move_to(temp_dot2))
            else:
                temp_group.add(green_glow.copy().move_to(temp_dot1), green_glow.copy().move_to(temp_dot2))
            temp_group.rotate(-22.5 * DEGREES, about_point = circle_center)
            return temp_group
        n = 50
        points_group = VGroup(*[get_two_points(circle_group[j]) for _ in range(n)])
        self.play(
            ShowSubmobjectsOneByOne(points_group),
            run_time = 10,
            rate_func = linear
        )
        self.wait()

        density_function = MathTex("f(\\theta)", "=", "\\frac{\\sin\\theta\\sinh^{-1}(\\tan\\theta)}{\\pi/2}")
        probability_statements = VGroup(
            MathTex("\\mathbb{P}(E)", "=", "\\int_0^{\\pi/2}\\mathbb{P}(E|\\theta)f(\\theta)\\,d\\theta"),
            MathTex("=", "\\frac{427+60\\pi^2-480\\ln 2}{180\\pi^2}"),
            MathTex("\\approx", "0.3864")
        )
        probability_statements[1].scale(1 - 1 / 8)
        density_function.move_to(circle_group[j].get_center() + 6 * UP / 4 + 6.75 * RIGHT)
        box = SurroundingRectangle(density_function)
        probability_statements[0].set_background_stroke(width = 5)
        probability_statements[0].move_to(circle_group[j].get_center() - 1 * UP / 4 + 6.75 * RIGHT)
        probability_statements[1].move_to(probability_statements[0][1].get_center() + probability_statements[1].get_center() - probability_statements[1][0].get_center())
        probability_statements[2].move_to(probability_statements[0][1].get_center() + probability_statements[2].get_center() - probability_statements[2][0].get_center())
        probability_statements[1].shift(3 * DOWN / 2)
        probability_statements[2].shift(6 * DOWN / 2)

        self.play(
            Write(density_function),
            Create(box),
            #event_definition.animate.shift(5 * RIGHT)
        )
        self.wait()
        self.play(Write(probability_statements[0]))
        self.wait()
        self.play(Write(probability_statements[1:]))
        self.wait()
        self.play(
            self.camera.frame.animate.move_to(ORIGIN).scale(2),
            FadeOut(VGroup(density_function, random_chord, points_group, theta_text, tangent_angle, tangent_line, method_names[j], box, probability_statements))
        )
        self.wait()

if __name__ == "__main__":
    module_name = os.path.abspath(__file__)
    #output_location = "C:\ManimCE\media"
    #clear_cmd = "cls"
    #command_A = "manim " + module_name + " " + "RealCase" + " " + "-pql -n 42" + " --media_dir " + output_location
    #command_A = "manim " + module_name + " --media_dir " + output_location + " -pqh"
    #command_A = "manim "+ "-p" + " " + module_name + " " + "SphericalTautochrone" + " -n 0"+ " --renderer=opengl"
    #command_A = "manim "+ "-pqh" + " " + module_name + " " + "RandomPointRandomAngle3" + " -n 13,17"
    #command_A = "manim "+ "-pqh" + " " + module_name + " " + "CircumferenceMethodProbability1" + " -n 0,5" + "  --disable_caching"
    command_A = "manim "+ "-pqh" + " " + module_name + " --disable_caching"
    #command_A = "manim "+ "-pqh" + " " + module_name
    #command_A = "manim "+ "--write_to_movie" + " " + module_name + " " + "SphericalTautochrone" + " -n 0"+ " --renderer=opengl"
    #os.system(clear_cmd)
    os.system(command_A)