; ================================
; COMMENTS
; ================================
(comment) @comment

((comment) @comment.documentation
 (#match? @comment "^/\\*\\*"))

; ================================
; STRUCTURE
; ================================
(CONFIG) @keyword.directive
(LUA) @keyword.directive
(EVENT) @keyword
(SECTION
  "SECTION" @keyword.directive
  name: (_) @label)

; ================================
; CONFIG
; ================================
(CONFIG
  value: (_) @number)

; ================================
; EVENTS (semantic types)
; ================================
(note_event) @type.builtin
(rest_event) @type.builtin
(ptech_event) @type.builtin
(utech_event) @type.builtin
(chord_event) @type.builtin
(trill_event) @type.builtin
(lua_event) @type.builtin

; ================================
; EVENT FIELDS
; ================================
(pitch) @string
[
  (ptech_event
    technique: (identifier) @string)

  (ptech_event
    techniques: (technique_group
      technique: (identifier) @string))

  (utech_event
    technique: (identifier) @string)

  (utech_event
    techniques: (technique_group
      technique: (identifier) @string))
]

[
  (note_event
    duration: (number) @variable.parameter)
  (rest_event
    duration: (number) @variable.parameter)
  (ptech_event
    duration: (number) @variable.parameter)
  (utech_event
    duration: (number) @variable.parameter)
]

; ================================
; ACTIONS (runtime layer)
; ================================
(exec receiver: (identifier) @type)
(exec args: (pdargs) @number)

; actions keywords
(exec "sendto" @function)
(exec "luacall" @function)
(delay "delay" @function)

; delay
;(delay amount: (number) @variable.parameter)
(delay amount: (number) @variable.parameter)
(delay unit: (time_unit) @variable.parameter)
