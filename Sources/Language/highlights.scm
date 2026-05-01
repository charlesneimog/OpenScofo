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
(EVENT) @keyword

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
  (utech_event
    technique: (identifier) @string)
  (ptech_event
    technique: (identifier) @string)
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
(exec receiver: (identifier) @variable.member)
(exec args: (pdargs) @number)

; actions keywords
(exec "sendto" @keyword)
(exec "luacall" @keyword)
(delay "delay" @keyword)

; delay 
(delay amount: (number) @variable.parameter)
(delay unit: (time_unit) @type)
